from pathlib import Path

def test_outputs():
    outfile1=Path("../src/plot_dsig_dtheta_vs_theta.png")
    if outfile1.is_file():
        print("file1 ok")
    outfile2=Path("../src/plot_theta_vs_b.png")
    if outfile2.is_file():
        print("file2 ok")
    outfile3=Path("../src/plot_b_vs_theta.png")
    if outfile3.is_file():
        print("file3 ok")
