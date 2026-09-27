import csv
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from illumina_qc.sav import instrument_type, load_summary_csv


class SAVTests(unittest.TestCase):
    def test_published_summary_families(self):
        records = list(load_summary_csv(ROOT / 'examples/output/illumina_runs_sav_summary.csv'))
        self.assertEqual(len(records), 292)
        self.assertEqual({name: sum(summary['Instrument'] == name for _, summary, _ in records) for name in ('MiSeq', 'NextSeq 550', 'NextSeq 1000/2000')}, {'MiSeq': 193, 'NextSeq 550': 79, 'NextSeq 1000/2000': 20})
        self.assertTrue(all(lots is None and not indexes for lots, _, indexes in records))
        self.assertEqual(instrument_type(records[0][1]['RunID']), 'MiSeq')

    def test_reload_is_idempotent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            command = [sys.executable, str(ROOT / 'src/illumina_qc/sav.py'), '--summary-csv', str(ROOT / 'examples/output/illumina_runs_sav_summary.csv'), '--database', str(root / 'metrics.sqlite'), '--export-csv', str(root / 'all.csv')]
            subprocess.run(command, check=True, capture_output=True)
            subprocess.run(command, check=True, capture_output=True)
            with sqlite3.connect(root / 'metrics.sqlite') as connection:
                self.assertEqual(connection.execute('SELECT COUNT(*) FROM summary').fetchone()[0], 292)
                self.assertEqual(connection.execute("SELECT name FROM sqlite_master WHERE name IN ('run_lots', 'percent_index')").fetchall(), [])
            for name, expected in (('miseq.csv', 193), ('nextseq_550.csv', 79), ('nextseq_1000_2000.csv', 20)):
                with (root / name).open(newline='') as handle:
                    self.assertEqual(sum(1 for _ in csv.DictReader(handle)), expected)


if __name__ == '__main__':
    unittest.main()
