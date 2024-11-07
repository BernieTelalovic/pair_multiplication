import unittest

if __name__ == '__main__':
    loader = unittest.TestLoader()
    suite = loader.discover('tests')  # Update path as needed
    runner = unittest.TextTestRunner()
    runner.run(suite)
