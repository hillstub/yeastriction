import os
from dotenv import load_dotenv

load_dotenv()

DATABASE_URL = os.getenv('DATABASE_URL', 'sqlite:///./data/database.db')
GENOMES_DIR = os.getenv('GENOMES_DIR', './data')
ALLOW_IMPORT = os.getenv('ALLOW_IMPORT', 'False').lower() == 'true'