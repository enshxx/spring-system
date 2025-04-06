from pathlib import Path
from csv_plt import show_plt


script_dir = Path(__file__).parent

csv_path = script_dir / "csv_data" / "integrate_data.csv"

show_plt(csv_path)
