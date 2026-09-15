#!/usr/bin/env python3
import evaluate_formal
from transport_ratio_model import load_bundle
from v4_common import evaluation_roles, load_rows

evaluate_formal.load_bundle = load_bundle
evaluate_formal.formal_evaluation_roles = evaluation_roles
evaluate_formal.load_rows = load_rows

if __name__ == "__main__":
    evaluate_formal.main()
