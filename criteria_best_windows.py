#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


def parse_percent(value: str) -> float:
    return float(value.strip().replace("%", ""))


def detect_section(row: list[str]) -> str | None:
    if len(row) > 1 and row[1].strip() == "Photoinduction":
        return "Photoinduction"
    if len(row) > 1 and row[1].strip() == "Pomeron":
        return "Pomeron"
    return None


def is_observable_row(row: list[str]) -> bool:
    return len(row) > 1 and row[1].strip() not in {"", "Observable"}


def load_ranking(input_path: Path) -> dict[str, dict[str, dict[str, str]]]:
    data: dict[str, dict[str, dict[str, str]]] = {"Photoinduction": {}, "Pomeron": {}}
    section = None

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            new_section = detect_section(row)
            if new_section is not None:
                section = new_section
                continue

            if section is None or not is_observable_row(row):
                continue

            observable = row[1].strip()
            data[section][observable] = {
                "eta_window": row[2].strip(),
                "energy": row[3].strip(),
                "ks_window_value": row[4].strip(),
                "ks_integrated_value": row[5].strip(),
                "improvement": row[6].strip(),
                "cross_ks_value": row[7].strip(),
                "cross_improvement": row[8].strip(),
            }

    return data


def choose_best_ranges(data: dict[str, dict[str, dict[str, str]]]) -> list[dict[str, str]]:
    photo = data["Photoinduction"]
    pomeron = data["Pomeron"]
    observables = sorted(set(photo) & set(pomeron))
    rows: list[dict[str, str]] = []

    for observable in observables:
        photo_row = photo[observable]
        pomeron_row = pomeron[observable]

        combined_if_photo = (
            parse_percent(photo_row["improvement"]) +
            parse_percent(pomeron_row["cross_improvement"])
        )
        combined_if_pomeron = (
            parse_percent(photo_row["cross_improvement"]) +
            parse_percent(pomeron_row["improvement"])
        )

        if combined_if_photo >= combined_if_pomeron:
            rows.append(
                {
                    "Observable": observable,
                    "Selected_from": "Photoinduction",
                    "Eta window": photo_row["eta_window"],
                    "Energy": photo_row["energy"],
                    "Photoinduction improvement %": f"{parse_percent(photo_row['improvement']):.2f}",
                    "Pomeron improvement %": f"{parse_percent(pomeron_row['cross_improvement']):.2f}",
                    "Combined improvement %": f"{combined_if_photo:.2f}",
                    "Selection criterion": "Maximum combined improvement across both tables",
                }
            )
        else:
            rows.append(
                {
                    "Observable": observable,
                    "Selected_from": "Pomeron",
                    "Eta window": pomeron_row["eta_window"],
                    "Energy": pomeron_row["energy"],
                    "Photoinduction improvement %": f"{parse_percent(photo_row['cross_improvement']):.2f}",
                    "Pomeron improvement %": f"{parse_percent(pomeron_row['improvement']):.2f}",
                    "Combined improvement %": f"{combined_if_pomeron:.2f}",
                    "Selection criterion": "Maximum combined improvement across both tables",
                }
            )

    return rows


def write_output(output_path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "Observable",
        "Selected_from",
        "Eta window",
        "Energy",
        "Photoinduction improvement %",
        "Pomeron improvement %",
        "Combined improvement %",
        "Selection criterion",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    default_input = Path("/home/guille/Downloads/Ranking_singletop.csv")
    default_output = Path("/home/guille/Downloads/Best_common_ranges_singletop.csv")

    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_input
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else default_output

    data = load_ranking(input_path)
    rows = choose_best_ranges(data)
    write_output(output_path, rows)

    print(f"Wrote {len(rows)} rows to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
