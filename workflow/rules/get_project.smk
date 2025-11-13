# Determine TCGA projects dynamically from outputs of download script,
# falling back to available directories or a static list.
def get_projects():
    csv_path = "data/processed/gdc_pancan/normal_counts_by_project.csv"
    if os.path.exists(csv_path):
        try:
            with open(csv_path, newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                projects = [row[0] for row in reader if row and row[0]]
                projects = sorted(set(projects))
                if projects:
                    return projects
        except Exception:
            pass

    # Fallback: infer from downloaded GDC directories
    dirs = [
        os.path.basename(p.rstrip("/"))
        for p in glob.glob("data/raw/GDCdata/TCGA-*")
        if os.path.isdir(p)
    ]
    dirs = sorted(set(dirs))
    if dirs:
        return dirs