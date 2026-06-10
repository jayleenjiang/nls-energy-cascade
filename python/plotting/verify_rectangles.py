import pandas as pd
from fractions import Fraction
import math 

def format_rational(x):
    if isinstance(x, str) and '/' in x:
        try:
            num, den = map(int, x.split('/'))
            return Fraction(num, den)
        except ValueError:
            return Fraction(0) 
    try:
        return Fraction(x)
    except (ValueError, TypeError):
        return Fraction(0)

def verify_file(filename):
    try:
        df = pd.read_csv(filename)
    except (FileNotFoundError, pd.errors.EmptyDataError) as e:
        print(f"Error reading file: {e}")
        return False

    all_rectangles = True
    num_families = df['family_id'].nunique()

    for family_id, group in df.groupby('family_id'):
        if len(group) != 4:
            print(f"Family {family_id}: Incorrect number of points ({len(group)}), skipping.")
            all_rectangles = False
            continue

        parents = group[group['point_type'] == 'parent']
        children = group[group['point_type'] == 'child']

        if len(parents) != 2 or len(children) != 2:
             print(f"Family {family_id}: Incorrect parent/child count, skipping.")
             all_rectangles = False
             continue

        P1_re = format_rational(parents.iloc[0]['real'])
        P1_im = format_rational(parents.iloc[0]['imag'])
        P2_re = format_rational(parents.iloc[1]['real'])
        P2_im = format_rational(parents.iloc[1]['imag'])
        C1_re = format_rational(children.iloc[0]['real'])
        C1_im = format_rational(children.iloc[0]['imag'])
        C2_re = format_rational(children.iloc[1]['real'])
        C2_im = format_rational(children.iloc[1]['imag'])

        MidP_re = (P1_re + P2_re) / 2
        MidP_im = (P1_im + P2_im) / 2
        MidC_re = (C1_re + C2_re) / 2
        MidC_im = (C1_im + C2_im) / 2

        midpoints_match = (MidP_re == MidC_re and MidP_im == MidC_im)

        DistP_sq = (P1_re - P2_re)**2 + (P1_im - P2_im)**2
        DistC_sq = (C1_re - C2_re)**2 + (C1_im - C2_im)**2

        lengths_match = (DistP_sq == DistC_sq)

        if not (midpoints_match and lengths_match):
            print(f"Family {family_id}")
            print(f"  MidP: ({float(MidP_re):.6f}, {float(MidP_im):.6f}i)")
            print(f"  MidC: ({float(MidC_re):.6f}, {float(MidC_im):.6f}i)")
            print(f"  DistP Sq: {float(DistP_sq):.6f}")
            print(f"  DistC Sq: {float(DistC_sq):.6f}")
            all_rectangles = False

    if all_rectangles:
        print(f"Success.")
    else:
        print(f"Failure.")

    return all_rectangles

if __name__ == '__main__':
    verify_file("rectangles_N5_R100_f1_to_f2.csv")
    verify_file("rectangles_N5_R100_f2_to_f3.csv")
    verify_file("rectangles_N5_R100_f3_to_f4.csv")
    verify_file("rectangles_N5_R100_f4_to_f5.csv")


