"""Reader for NLS3_gedmd_stream output (OUT.bin + OUT.json)."""
import json
import numpy as np


def read(prefix):
    meta = json.load(open(prefix + ".json"))
    K, G, NL = meta["K"], meta["groups"], meta["nlag"]
    nobs = len(meta["observables"])
    groups = []
    with open(prefix + ".bin", "rb") as f:
        assert f.read(8) == b"NLS3GED1"
        hdr = np.frombuffer(f.read(24), "<i8")
        assert tuple(hdr) == (K, G, NL), (hdr, K, G, NL)
        for _ in range(G):
            n = np.frombuffer(f.read(8), "<f8")[0]
            A = np.frombuffer(f.read(8 * K * K), "<f8").reshape(K, K, order="F")
            B = np.frombuffer(f.read(8 * K * K), "<f8").reshape(K, K, order="F")
            sL = np.frombuffer(f.read(8 * K), "<f8")
            sL2 = np.frombuffer(f.read(8 * K), "<f8")
            sx = np.frombuffer(f.read(40), "<f8")
            so = np.frombuffer(f.read(8 * nobs), "<f8")
            acf = np.frombuffer(f.read(8 * nobs * NL), "<f8").reshape(nobs, NL)
            acfn = np.frombuffer(f.read(8 * NL), "<f8")
            groups.append(dict(n=n, A=A, B=B, sL=sL, sL2=sL2, sx=sx, so=so, acf=acf, acfn=acfn))
    terms = [((t[0], t[1], t[2]), (t[3], t[4], t[5]), t[6], t[7], "s" if t[8] else "c") for t in meta["terms"]]
    return meta, terms, groups


def combine(groups, which):
    n = sum(groups[i]["n"] for i in which)
    out = {"n": n}
    for key in ("A", "B", "sL", "sL2", "sx", "so", "acf", "acfn"):
        out[key] = sum(groups[i][key] for i in which)
    out["A"] = out["A"] / n; out["B"] = out["B"] / n
    out["mL"] = out["sL"] / n; out["sdL"] = np.sqrt(np.maximum(out["sL2"] / n - out["mL"] ** 2, 0))
    out["mean_x"] = out["sx"] / n
    return out


def sub_index(terms, d, M):
    """Indices of the nested sub-dictionary (d, M), in the same order make_dictionary(d, M) would produce."""
    import gedmd_n3 as G
    pos = {t: k for k, t in enumerate(terms)}
    return np.array([pos[t] for t in G.make_dictionary(d, M)])
