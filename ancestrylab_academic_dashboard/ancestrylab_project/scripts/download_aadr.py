#!/usr/bin/env python3
"""Metadata-first, resumable and checksum-verified AADR Dataverse downloader."""
import argparse, hashlib, json, os, pathlib, shutil, sys, time
import requests

PID = "doi:10.7910/DVN/FFIDCW"
BASE = "https://dataverse.harvard.edu"
S = requests.Session(); S.headers.update({"User-Agent":"AncestryLab/0.1 (reproducible research)", "Accept":"application/json"})

def dataset():
    r = S.get(f"{BASE}/api/datasets/:persistentId/", params={"persistentId": PID}, timeout=(20,60)); r.raise_for_status()
    obj = r.json()
    if obj.get("status") != "OK": raise RuntimeError(f"Dataverse API returned {obj.get('status')}")
    return obj["data"]["latestVersion"]

def select_files(version, panel):
    all_files = [x["dataFile"] for x in version.get("files", [])]
    selected = {}
    for ext in ("ind", "snp", "geno", "anno"):
        matches = [f for f in all_files if f.get("filename", "").lower().endswith("."+ext)
                   and panel.lower() in f.get("filename", "").lower()
                   and "pub" in f.get("filename", "").lower()]
        if ext in ("ind", "snp", "geno"):
            patched = [f for f in matches if "patch" in f.get("filename", "").lower()]
            if patched: matches = patched
        if not matches: raise FileNotFoundError(f"No public AADR {panel} .{ext} file in the latest version")
        matches.sort(key=lambda x: int(x.get("filesize",0)), reverse=True)
        if len(matches)>1 and matches[0].get("filesize")==matches[1].get("filesize"):
            raise RuntimeError(f"Ambiguous {ext} candidates: {[m['filename'] for m in matches]}")
        selected[ext]=matches[0]
    return selected

def md5(path):
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def download(info, folder):
    folder.mkdir(parents=True, exist_ok=True); dest=folder/info["filename"]; part=folder/(info["filename"]+".part")
    expected=int(info.get("filesize") or 0); expected_md5=(info.get("md5") or "").lower()
    if dest.exists() and (not expected or dest.stat().st_size==expected):
        if expected_md5 and md5(dest)!=expected_md5: dest.unlink()
        else: print(f"Verified existing: {dest.name}"); return dest
    free=shutil.disk_usage(folder).free
    if expected and free < expected + 512*1024*1024: raise OSError(f"Insufficient disk: need {expected/1e9:.2f} GB plus 0.5 GB reserve")
    offset=part.stat().st_size if part.exists() else 0
    headers={"Range":f"bytes={offset}-"} if offset else {}
    with S.get(f"{BASE}/api/access/datafile/{info['id']}", headers=headers, stream=True, timeout=(30,180)) as r:
        if offset and r.status_code==206: mode="ab"
        else: offset=0; mode="wb"; r.raise_for_status()
        total=expected or None; done=offset; last=time.time()
        with open(part,mode) as out:
            for chunk in r.iter_content(8*1024*1024):
                if chunk:
                    out.write(chunk); done+=len(chunk)
                    if time.time()-last>3:
                        print(f"{dest.name}: {done/1e9:.2f}/{(total or done)/1e9:.2f} GB", flush=True); last=time.time()
    if expected and part.stat().st_size!=expected: raise IOError(f"Incomplete download ({part.stat().st_size}/{expected}); rerun to resume")
    if expected_md5 and md5(part)!=expected_md5: raise IOError(f"MD5 mismatch for {dest.name}; partial file retained for inspection")
    part.replace(dest); print(f"Downloaded and verified: {dest}"); return dest

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--panel",choices=["1240K","HO","2M"],default="1240K"); ap.add_argument("--out",default="data/aadr"); ap.add_argument("--metadata-only",action="store_true")
    a=ap.parse_args(); v=dataset(); selected=select_files(v,a.panel); out=pathlib.Path(a.out)
    manifest={"dataset_pid":PID,"version_number":v.get("versionNumber"),"version_date":v.get("versionDate"),"panel":a.panel,"files":{k:{"id":x["id"],"filename":x["filename"],"size":x.get("filesize"),"md5":x.get("md5")} for k,x in selected.items()}}
    out.mkdir(parents=True,exist_ok=True); (out/"dataverse_manifest.json").write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))
    if a.metadata_only:
        # Small files only: enough for exact sample/population discovery.
        for key in ("ind", "anno"): download(selected[key], out)
        from collections import Counter
        counts = Counter()
        with open(out / selected["ind"]["filename"], encoding="utf-8") as fh:
            for line in fh:
                fields = line.split()
                if len(fields) >= 3: counts[fields[2]] += 1
        import csv
        with open(out / "population_inventory.csv", "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh); writer.writerow(["population", "n_individuals"])
            writer.writerows(sorted(counts.items(), key=lambda x: x[0]))
        print(f"Wrote exact population inventory: {out / 'population_inventory.csv'}")
    else:
        for key in ("ind","snp","anno","geno"): download(selected[key],out)
if __name__=="__main__": main()
