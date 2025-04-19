import unittest
import os
import logging
from Bio import SeqIO
from filter import filter_fastq

class TestFastqFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_input = "test_input.fastq"
        cls.test_output = "test_output.fastq"
        cls.test_log = "test.log"
        
        # Создаем FASTQ для теста
        with open(cls.test_input, 'w') as f:
            f.write("@seq1\nACGT\n+\nIIII\n"         # 50% GC, длина 4, качество ~40
                    "@seq2\nAAAA\n+\nIIII\n"         # 0% GC, длина 4, качество ~40
                    "@seq3\nTGCA\n+\n!!!!\n"         # 50% GC, длина 4, качество ~0
                    "@seq4\nGGGGGG\n+\nIIIIII\n")    # 100% GC, длина 6, качество ~40

    def setUp(self):
        # Очищаем выходные файлы
        for f in [self.test_output, self.test_log]:
            if os.path.exists(f):
                os.remove(f)
        
        # Настраиваем логгер
        logging.basicConfig(
            level=logging.INFO,
            filename=self.test_log,
            filemode='w',
            force=True
        )

    def test_basic_filtering(self):
        """Без фильтров должны пройти все записи"""
        count = filter_fastq(self.test_input, self.test_output)
        self.assertEqual(count, 4)
        with open(self.test_output) as f:
            records = list(SeqIO.parse(f, "fastq"))
            self.assertEqual(len(records), 4)

    def test_gc_filter_lower(self):
        """Фильтрация по минимальному GC=50%"""
        count = filter_fastq(self.test_input, self.test_output, gc_bounds=(50, 100))
        self.assertEqual(count, 3)  # ACGT, TGCA, GGGGGG
        with open(self.test_output) as f:
            records = list(SeqIO.parse(f, "fastq"))
            seqs = {str(r.seq) for r in records}
            self.assertEqual(seqs, {"ACGT", "TGCA", "GGGGGG"})

    def test_gc_filter_upper(self):
        """Фильтрация по максимальному GC=0%"""
        count = filter_fastq(self.test_input, self.test_output, gc_bounds=0)
        self.assertEqual(count, 1)  # только AAAA
        with open(self.test_output) as f:
            records = list(SeqIO.parse(f, "fastq"))
            self.assertEqual(str(records[0].seq), "AAAA")

    def test_quality_filter(self):
        """Фильтрация по качеству >=30"""
        count = filter_fastq(self.test_input, self.test_output, quality_threshold=30)
        self.assertEqual(count, 3)  # ACGT, AAAA, GGGGGG (TGCA имеет качество ~0)
        with open(self.test_output) as f:
            records = list(SeqIO.parse(f, "fastq"))
            seqs = {str(r.seq) for r in records}
            self.assertEqual(seqs, {"ACGT", "AAAA", "GGGGGG"})

    def test_length_filter(self):
        """Фильтрация по длине >=5"""
        count = filter_fastq(self.test_input, self.test_output, length_bounds=(5, 100))
        self.assertEqual(count, 1)  # только GGGGGG (длина 6)
        with open(self.test_output) as f:
            records = list(SeqIO.parse(f, "fastq"))
            self.assertEqual(len(records), 1)
            self.assertEqual(str(records[0].seq), "GGGGGG")

    def test_file_creation(self):
        """Проверка создания выходного файла"""
        self.assertFalse(os.path.exists(self.test_output))
        filter_fastq(self.test_input, self.test_output)
        self.assertTrue(os.path.exists(self.test_output))

    def test_logging(self):
        """Проверка логирования"""
        filter_fastq(self.test_input, self.test_output)
        self.assertTrue(os.path.exists(self.test_log))
        with open(self.test_log) as f:
            log = f.read()
            self.assertIn("Starting filtering", log)
            self.assertIn("Filtered 4 records", log)

    def test_error_handling(self):
        """Проверка обнаружения ошибки, если нет инпута"""
        with self.assertRaises(FileNotFoundError):
            filter_fastq("nonexistent.fastq", self.test_output)

    @classmethod
    def tearDownClass(cls):
        # Очистка тестовых файлов
        for f in [cls.test_input, cls.test_output, cls.test_log]:
            if os.path.exists(f):
                os.remove(f)

if __name__ == '__main__':
    unittest.main()

# pytest test_fastq_filter.py -v