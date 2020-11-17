"""
Written by Erin Molloy (molloy.erin.k@gmail.com) in Summer 2019.
"""
import argparse
import numpy
import os


def set_indelible_params(ngens, freqs, rates, alpha, sqlen, outf):
    """
    Parameters
    ----------
    ngens : int
            Number of genes
    freqs : list of floats (length 4)
            Alpha values (A, C, G, T) to define a Dirichlet distribution from
            which GTR base frequencies will be drawn.
    rates : list of floats (length 6)
            Alpha values (AC, AG, AT, CG, CT, GT) to define a Dirichlet
            distribution from which GTR transition rates will be drawn.
    alpha : list of floats (length 2)
            Meanlog and sdlog to define a Lognormal distribution from which
            alpha (to define a gamma for site rate heterogeneity) will be drawn.
    sqlen : int
            Sequence length
    outf : string
           output file name

    Returns
    -------
    Nothing, writes an output file
    """
    with open(outf, 'w') as f:
        f.write("GENE,SQLN,ALPH,fA,fC,fG,fT,rAC,rAG,rAT,rCG,rCT,rGT\n")
        fmat = numpy.random.dirichlet((freqs[0], freqs[1], freqs[2], freqs[3]),
                                      ngens)
        rmat = numpy.random.dirichlet((rates[0], rates[1], rates[2], rates[3],
                                       rates[4], rates[5]), ngens)
        avec = numpy.random.lognormal(alpha[0], alpha[1], ngens)

        gene = 1
        for j in range(ngens):
            fj = fmat[j, :]
            rj = rmat[j, :]
            f.write("%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                    % (gene, sqlen, avec[j],
                       fj[0], fj[1], fj[2], fj[3],
                       rj[0], rj[1], rj[2], rj[3], rj[4], rj[5]))
            gene = gene + 1


def main(args):
    set_indelible_params(args.ngens, args.freqs, args.rates, args.alpha,
                         args.sqlen, args.output)

    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--ngens", type=int, required=True,
                        help="Specify number of genes.")
    parser.add_argument("-f", "--freqs", type=float, nargs=4, required=True,
                        help="GTR base frequencies: Specify alpha parameters "
                             "(A, C, G, T) for Dirichlet distribution.")
    parser.add_argument("-r", "--rates", type=float, nargs=6, required=True,
                        help="GTR transition rate matrix: Specify alpha "
                             "parameters (AC, AG, AT, CG, CT, GT) for Dirichlet"
                             " distribution.")
    parser.add_argument("-a", "--alpha", type=float, nargs=2, required=True,
                        help="Alpha: Specify meanlog and sdlog for Lognormal "
                             "distribution.")
    parser.add_argument("-l", "--sqlen", type=int, required=True,
                        help="Specify sequence length (constant)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file")

    main(parser.parse_args())
