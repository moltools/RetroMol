import argparse
import json
import random
from tqdm import tqdm
from typing import Any

import aiohttp
import asyncio

from retromol.io.json import iter_json
from retromol.model.result import Result


# NPClassifier API settings
NP_API_BASES = [
    "https://npclassifier.ucsd.edu/classify",
    "https://npclassifier.gnps2.org/classify",
]

NP_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "BioNexus/0.2 (NPClassifier bulk annotate)",
}


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True)
    parser.add_argument("--out-tsv", required=True)

    parser.add_argument("--concurrency", type=int, default=4)
    parser.add_argument("--rate-per-sec", type=float, default=1.0)
    parser.add_argument("--timeout-s", type=int, default=30)
    parser.add_argument("--retries", type=int, default=3)

    return parser.parse_args()


def _as_text(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, list):
        return ";".join(map(str, x))
    return str(x)


async def _np_fetch_one(
        session: aiohttp.ClientSession,
        smiles: str,
        *,
        timeout_s: int,
        retries: int,
        delay_between: float | None,
        semaphore: asyncio.Semaphore,
) -> dict[str, Any] | None:
    """
    Fetch NPClassifier annotation for a single SMILES string.

    :param session: aiohttp client session
    :param smiles: SMILES string to query
    :param timeout_s: request timeout in seconds
    :param retries: number of retries on failure
    :param delay_between: delay between requests in seconds
    :param semaphore: asyncio semaphore for limiting concurrency
    :return: NPClassifier annotation blob or None
    """
    async with semaphore:
        if delay_between:
            await asyncio.sleep(delay_between)

        backoff = 0.5
        last_err: Exception | None = None

        for attempt in range(retries + 1):
            for base in NP_API_BASES:
                try:
                    async with session.get(
                            base,
                            params={"smiles": smiles},
                            timeout=timeout_s,
                            headers=NP_HEADERS,
                    ) as resp:
                        status = resp.status
                        text = await resp.text()
                        if status == 200:
                            # Header sometimes wrong; try parse anyway
                            try:
                                return json.loads(text)
                            except json.JSONDecodeError:
                                pass

                        if status == 429:
                            ra = resp.headers.get("Retry-After")
                            wait = float(ra) if (ra and ra.isdigit()) else backoff * (1 + random.random())
                            await asyncio.sleep(wait)
                            continue

                        if 500 <= status < 600:
                            last_err = RuntimeError(f"{base} returned {status}")
                            continue

                        if 400 <= status < 500:
                            print(f"NPClassifier 4xx ({status}) for SMILES={smiles}", status, smiles)
                            return None

                except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                    # Log and retry
                    last_err = e

            if attempt < retries:
                # Exponential backoff with jitter
                await asyncio.sleep(backoff * (1 + random.random()))
                backoff *= 2

        if last_err:
            print(f"NPClassifier giving up for SMILES={smiles} :: {last_err}")

        return None


async def np_fetch_many(
    smiles_list: list[str],
    *,
    concurrency: int,
    rate_per_sec: float | None,
    timeout_s: int,
    retries: int,
    pbar: tqdm | None = None,
) -> dict[str, dict[str, Any] | None]:
    """
    Fetch NPClassifier annotations for multiple SMILES strings concurrently.

    :param smiles_list: list of SMILES strings to query
    :param concurrency: maximum number of concurrent requests
    :param rate_per_sec: maximum request rate per second
    :param timeout_s: request timeout in seconds
    :param retries: number of retries on failure
    :param pbar: optional tqdm progress bar
    :return: dictionary mapping SMILES strings to their NPClassifier annotation blobs or None
    """
    semaphore = asyncio.Semaphore(max(1, concurrency))
    delay_between = (1.0 / rate_per_sec) if rate_per_sec and rate_per_sec > 0 else None

    connector = aiohttp.TCPConnector(
        limit=max(4, concurrency),
        limit_per_host=max(2, concurrency // 2),
        enable_cleanup_closed=True,
    )
    timeout = aiohttp.ClientTimeout(total=None)

    out: dict[str, dict[str, Any] | None] = {}
    async with aiohttp.ClientSession(timeout=timeout, connector=connector, headers=NP_HEADERS) as session:
        async def _run(smi: str):
            out[smi] = await _np_fetch_one(
                session, smi,
                timeout_s=timeout_s,
                retries=retries,
                delay_between=delay_between,
                semaphore=semaphore,
            )
            if pbar is not None:
                pbar.update(1)
        tasks = [asyncio.create_task(_run(smi)) for smi in smiles_list]
        await asyncio.gather(*tasks)
    return out




def main() -> None:
    args = cli()

    rows = []
    smiles_list = []

    for d_idx, d in tqdm(enumerate(iter_json(args.jsonl, jsonl=True)), desc="Reading results"):
        r = Result.from_dict(d)

        smiles = r.submission.smiles
        coverage = r.calculate_coverage()

        rows.append({
            "index": d_idx,
            "smiles": smiles,
            "coverage": coverage,
        })

        if smiles:
            smiles_list.append(smiles)

    unique_smiles = sorted(set(smiles_list))

    with tqdm(total=len(unique_smiles), desc="NPClassifier") as pbar:
        annotations = asyncio.run(
            np_fetch_many(
                unique_smiles,
                concurrency=args.concurrency,
                rate_per_sec=args.rate_per_sec,
                timeout_s=args.timeout_s,
                retries=args.retries,
                pbar=pbar,
            )
        )

    with open(args.out_tsv, "w") as fout:
        fout.write(
            "index\tsmiles\tcoverage\t"
            "np_pathway\tnp_superclass\tnp_class\tisglycoside\traw_json\n"
        )

        for row in rows:
            ann = annotations.get(row["smiles"])

            if ann is None:
                pathway = superclass = np_class = isglycoside = raw_json = ""
            else:
                pathway = _as_text(ann.get("pathway_results"))
                superclass = _as_text(ann.get("superclass_results"))
                np_class = _as_text(ann.get("class_results"))
                isglycoside = _as_text(ann.get("isglycoside"))
                raw_json = json.dumps(ann, sort_keys=True)

            fout.write(
                f"{row['index']}\t"
                f"{row['smiles']}\t"
                f"{row['coverage']}\t"
                f"{pathway}\t"
                f"{superclass}\t"
                f"{np_class}\t"
                f"{isglycoside}\t"
                f"{raw_json}\n"
            )



if __name__ == "__main__":
    main()
