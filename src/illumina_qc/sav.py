"""SAV ingestion: run XML, InterOp summary, SQLite tables, and Tableau exports."""

import argparse
import csv
import re
import sqlite3
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path


READ_FIELDS = {
    "Percent_Error_Rate": "error_rate",
    "Intensity_Cycle_1": "first_cycle_intensity",
    "Percent_Q30": "percent_gt_q30",
    "Legacy_Phasing_Rate": "phasing",
    "Legacy_Prephasing_Rate": "prephasing",
}
# These method names and labels come from the supplied sav_variables.py.
COL_READS = (("Total_Reads", "total_reads"), ("Total_Reads_PF", "total_pf_reads"),
             ("Percent_Reads_PF", "total_fraction_mapped_reads"),
             ("CV_Reads", "mapped_reads_cv"), ("Min_Reads", "min_mapped_reads"),
             ("Max_Reads", "max_mapped_reads"))
COL_SUM1 = (("Yield_Total", "yield_g"), ("Percent_Aligned", "percent_aligned"),
            ("Percent_Occupied_Wells", "percent_occupied"))
COL_SUM2 = (("Density", "density"), ("PF_Percent", "percent_pf"))
INDEX_FIELDS = {
    "Index_Number": "id", "Sample_ID": "sample_id", "Project": "project_name",
    "Index_1_I7": "index1", "Index_2_I5": "index2",
    "Percent_Reads_Identified": "fraction_mapped",
}


def xml_value(root, *names):
    wanted = {name.casefold() for name in names}
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1].casefold() in wanted and element.text and element.text.strip():
            return element.text.strip()
    return None


def xml_section(root, name):
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1].casefold() == name.casefold():
            return element
    return None


def run_date(folder_name):
    match = re.match(r"^(\d{6})_", folder_name)
    return datetime.strptime(match.group(1), "%y%m%d").date().isoformat() if match else None


def instrument_type(run_id):
    parts = run_id.split("_")
    instrument = parts[1].upper() if len(parts) > 1 else ""
    if instrument.startswith("VL"):
        return "NextSeq 1000/2000"
    if instrument.startswith("M"):
        return "MiSeq"
    if instrument.startswith(("NDX", "NS", "NB")):
        return "NextSeq 550"
    return None


def find_run_parameters(folder):
    return next((path for path in folder.iterdir() if path.name.casefold() == "runparameters.xml"), None)


def panel_from_samplesheet(folder):
    for path in folder.iterdir():
        if "samplesheet" in path.name.casefold() and path.suffix.casefold() == ".csv":
            content = path.read_text(encoding="utf-8-sig", errors="replace").casefold()
            for word, panel in (("myeloid", "Myeloid"), ("cftr", "CF"), ("exome", "WES")):
                if word in content:
                    return panel
    return None


def parse_run_parameters(folder):
    parameters = find_run_parameters(folder)
    if parameters is None:
        raise FileNotFoundError(f"No RunParameters.xml in {folder}")
    root = ET.parse(parameters).getroot()
    row = {
        "RunID": xml_value(root, "runid", "runidentifier") or folder.name,
        "RunName": xml_value(root, "experimentname", "runname") or "NA",
        "Panel": panel_from_samplesheet(folder),
        "Instrument": instrument_type(folder.name),
        "RunDate": run_date(folder.name),
        "RunNumber": xml_value(root, "runnumber"),
        "RTAVersion": xml_value(root, "rtaversion"),
        "Instrument_SN": xml_value(root, "instrumentserialnumber", "instrumentid", "scannerid"),
    }
    for tag, prefix in (("flowcellrfidtag", "FLOWCELL"), ("pr2bottlerfidtag", "PR2BOTTLE"), ("reagentkitrfidtag", "REAGENTKIT")):
        section = xml_section(root, tag)
        row[f"{prefix}_SerialNo"] = xml_value(section, "serialnumber", "serialno", "serial") if section is not None else None
        row[f"{prefix}_PartNo"] = xml_value(section, "partnumber", "partno", "partnumberstring") if section is not None else None
        row[f"{prefix}_Expiry"] = xml_value(section, "expirationdate", "expirydate", "expiration") if section is not None else None
    setup = xml_section(root, "setup")
    if setup is None:
        setup = root
    row["Software"] = xml_value(setup, "applicationname", "software", "softwarename")
    row["Software_Ver"] = xml_value(setup, "applicationversion", "softwareversion")
    return row


def run_params(folder):
    """Extract the original SAV run_lots fields from RunParameters XML."""
    return parse_run_parameters(folder)


def call(obj, name, default=None):
    method = getattr(obj, name, None)
    return method() if callable(method) else default


def mean(value):
    return value.mean() if hasattr(value, "mean") else value


def read_interop(folder):
    try:
        from interop import py_interop_run_metrics, py_interop_summary
    except ImportError as exc:
        raise RuntimeError("Instrument ingestion requires the optional Illumina interop Python package") from exc

    metrics = py_interop_run_metrics.run_metrics()
    metrics.read(str(folder))
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(metrics, summary)
    read_count = summary.size()
    reads = [summary.at(i).summary() for i in range(read_count)]
    total = summary.total_summary()
    index_lane = py_interop_summary.index_lane_summary()
    py_interop_summary.summarize_index_metrics(metrics, 0, index_lane)

    def lane_values(name):
        values = [mean(call(summary.at(i).at(lane), name)) for i in range(read_count) for lane in range(summary.lane_count())]
        values = [value for value in values if isinstance(value, (int, float))]
        return sum(values) / len(values) if values else None

    lane_metrics = {label: lane_values(method) for label, method in COL_SUM2}
    row = {
        **{label: call(index_lane, method) for label, method in COL_READS},
        **{label: call(total, method) for label, method in COL_SUM1},
        **lane_metrics,
    }
    density = lane_metrics["Density"]
    row["Density"] = density / 1000 if density is not None else None
    for output, source in READ_FIELDS.items():
        for i in range(1, 5):
            if output in ("Legacy_Phasing_Rate", "Legacy_Prephasing_Rate") and i not in (1, 4):
                continue
            row[f"{output}_R{i}"] = call(reads[i - 1], source) if i <= read_count else None
        if output not in ("Legacy_Phasing_Rate", "Legacy_Prephasing_Rate"):
            row[f"{output}_Total"] = call(total, source)

    indexes = py_interop_summary.index_flowcell_summary()
    py_interop_summary.summarize_index_metrics(metrics, indexes)
    sample_rows = []
    if indexes.size():
        lane = indexes.at(0)
        for i in range(lane.size()):
            entry = lane.at(i)
            sample_rows.append({column: call(entry, getter) for column, getter in INDEX_FIELDS.items()})
    return row, sample_rows


def run_summary(folder):
    """Extract run level metrics using the original SAV InterOp methods."""
    return read_interop(folder)[0]


def index_per_sample(folder):
    """Extract sample index rows using the original SAV InterOp methods."""
    return read_interop(folder)[1]


def rows_from_folder(folder):
    lots = run_params(folder)
    summary, indexes = read_interop(folder)
    metadata = {key: lots[key] for key in ("RunID", "RunName", "Panel", "Instrument", "RunDate")}
    return lots, metadata | summary, [metadata | index for index in indexes]


def load_summary_csv(path):
    """Reload the supplied real SAV summary snapshot without XML/InterOp files."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"RunID", "Instrument", "RunDate", "Density", "PF_Percent"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Summary CSV missing fields: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            summary = {key: None if value in ("", "NULL", "NA") else value for key, value in row.items()}
            for key, value in summary.items():
                if key not in {"RunID", "RunName", "Panel", "Instrument", "RunDate"} and value is not None:
                    try:
                        summary[key] = float(value)
                    except ValueError:
                        pass
            if summary["Instrument"] == "NextSeq":
                summary["Instrument"] = "NextSeq 550"
            if summary["Instrument"] not in {"MiSeq", "NextSeq 550", "NextSeq 1000/2000"}:
                raise ValueError(f"Unexpected instrument: {summary['Instrument']}")
            summary["RunDate"] = summary["RunDate"].split()[0]
            yield None, summary, []


def replace_run(connection, lots, summary, indexes):
    run_id = summary["RunID"]
    if not run_id or (lots is not None and run_id != lots["RunID"]) or any(row["RunID"] != run_id for row in indexes):
        raise ValueError("Every table must use the same nonempty RunID")
    for table, rows in (("run_lots", [lots] if lots is not None else []), ("summary", [summary]), ("percent_index", indexes)):
        columns = sorted({column for row in rows for column in row})
        if not columns:
            continue
        connection.execute(f'CREATE TABLE IF NOT EXISTS "{table}" ("RunID" TEXT)')
        existing = {row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')}
        for column in columns:
            if column not in existing:
                if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", column):
                    raise ValueError(f"Invalid column: {column}")
                connection.execute(f'ALTER TABLE "{table}" ADD COLUMN "{column}" {"REAL" if isinstance(next((row[column] for row in rows if row.get(column) is not None), None), (int, float)) else "TEXT"}')
        connection.execute(f'DELETE FROM "{table}" WHERE "RunID" = ?', (run_id,))
        names = ", ".join(f'"{column}"' for column in columns)
        placeholders = ", ".join("?" for _ in columns)
        connection.executemany(f'INSERT INTO "{table}" ({names}) VALUES ({placeholders})', [[row.get(column) for column in columns] for row in rows])


def export_csv(connection, path, instrument=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    query = 'SELECT * FROM "summary"'
    params = []
    if instrument:
        query += ' WHERE "Instrument" = ?'
        params.append(instrument)
    cursor = connection.execute(query + ' ORDER BY "RunDate", "RunID"', params)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([column[0] for column in cursor.description])
        writer.writerows(cursor)


def sav_metrics_to_db(connection, records):
    """Write the three original SAV tables, keyed by RunID."""
    count = 0
    for lots, summary, indexes in records:
        with connection:
            replace_run(connection, lots, summary, indexes)
        count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--runs-dir", type=Path, help="Directory of local MiSeq/NextSeq run folders (requires InterOp)")
    source.add_argument("--summary-csv", type=Path, help="Saved SAV summary CSV for MiSeq and NextSeq")
    parser.add_argument("--database", type=Path, default=Path("examples/output/run_metrics.sqlite"))
    parser.add_argument("--export-csv", type=Path, default=Path("examples/output/tableau_run_metrics.csv"))
    args = parser.parse_args()
    if args.summary_csv:
        records = load_summary_csv(args.summary_csv)
    else:
        records = (rows_from_folder(folder) for folder in sorted(args.runs_dir.iterdir()) if folder.is_dir() and instrument_type(folder.name) and find_run_parameters(folder))
    args.database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(args.database) as connection:
        count = sav_metrics_to_db(connection, records)
        if not count:
            raise SystemExit("No eligible runs found; CSV was not exported")
        export_csv(connection, args.export_csv)
        for instrument, filename in (("MiSeq", "miseq.csv"), ("NextSeq 550", "nextseq_550.csv"), ("NextSeq 1000/2000", "nextseq_1000_2000.csv")):
            export_csv(connection, args.export_csv.parent / filename, instrument)
    print(f"Loaded {count} SAV runs to {args.database}; exported Tableau CSVs in {args.export_csv.parent}")


if __name__ == "__main__":
    main()
