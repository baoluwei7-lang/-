import hashlib
from pathlib import Path
from datetime import datetime
import csv


def calculate_sha256(file_path, chunk_size=1024 * 1024):
    """
    分块计算文件 SHA256
    """
    sha256_hash = hashlib.sha256()

    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            sha256_hash.update(chunk)

    return sha256_hash.hexdigest()


if __name__ == "__main__":

    # 当前脚本路径 notebooks/
    script_dir = Path(__file__).resolve().parent

    # 项目根目录
    project_root = script_dir.parent

    # 原始数据目录
    data_dir = project_root / "data_raw"

    # 输出文件
    output_file = project_root / "data_manifest.tsv"


    if not data_dir.exists():
        raise FileNotFoundError(
            f"数据目录不存在:\n{data_dir}"
        )


    print(f"扫描目录: {data_dir}")

    records = []

    # 当前日期
    download_date = datetime.now().strftime("%Y-%m-%d")


    # 遍历 data_raw 所有文件
    for file_path in sorted(data_dir.iterdir()):

        # 跳过文件夹
        if not file_path.is_file():
            continue

        print(f"\n正在处理: {file_path.name}")

        sha256 = calculate_sha256(file_path)

        record = {
            "dataset": file_path.stem,
            "source_url": "-",
            "accession": "-",
            "download_date": download_date,
            "file_name": file_path.name,
            "size_bytes": file_path.stat().st_size,
            "sha256": sha256,
            "license": "-",
            "用途": "-"
        }

        records.append(record)


    # 写入 TSV
    with open(
        output_file,
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset",
                "source_url",
                "accession",
                "download_date",
                "file_name",
                "size_bytes",
                "sha256",
                "license",
                "用途"
            ],
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(records)


    print("\n================================")
    print("完成!")
    print(f"共记录 {len(records)} 个文件")
    print(f"manifest 输出:")
    print(output_file)
    print("================================")