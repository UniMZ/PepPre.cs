namespace UniMZ.PepPre;

// copied from CIAAW on 2024-11-03 (http://ciaaw.org)
public readonly record struct DefaultIsotopeTable
{
    public static readonly Isotope H = new(
        [1.0078250322, 2.0141017781],
        [0.999855, 0.000145]
    );

    public static readonly Isotope C = new(
        [12.0, 13.003354835],
        [0.9894, 0.0106]
    );

    public static readonly Isotope N = new(
        [14.003074004, 15.000108899],
        [0.996205, 0.003795]
    );

    public static readonly Isotope O = new(
        [15.994914619, 16.999131757, 17.999159613],
        [0.9975715, 0.0003835, 0.002045]
    );

    public static readonly Isotope S = new(
        [31.972071174, 32.97145891, 33.9678670, 0.0, 35.967081],
        [0.9485588, 0.0076305, 0.0436527, 0.0, 0.000158]
    );

    public static readonly Isotope[] HCNOS = [H, C, N, O, S];
}

// abundances copied from Proteome-pI 2.0 on 2025-07-23 (https://isoelectricpointdb2.mimuw.edu.pl/statistics.html#aa_stats, All)
// formula: H C N O S
public readonly record struct DefaultAminoAcidTable
{
    public static readonly AminoAcid A = new("Ala", [5, 3, 1, 1, 0], 0.0872);
    public static readonly AminoAcid C = new("Cys", [5, 3, 1, 1, 1], 0.0146);
    public static readonly AminoAcid D = new("Asp", [5, 4, 1, 3, 0], 0.0549);
    public static readonly AminoAcid E = new("Glu", [7, 5, 1, 3, 0], 0.0636);
    public static readonly AminoAcid F = new("Phe", [9, 9, 1, 1, 0], 0.0378);
    public static readonly AminoAcid G = new("Gly", [3, 2, 1, 1, 0], 0.0704);
    public static readonly AminoAcid H = new("His", [7, 6, 3, 1, 0], 0.0232);
    public static readonly AminoAcid I = new("Ile", [11, 6, 1, 1, 0], 0.0519);
    public static readonly AminoAcid K = new("Lys", [12, 6, 2, 1, 0], 0.0506);
    public static readonly AminoAcid L = new("Leu", [11, 6, 1, 1, 0], 0.0967);
    public static readonly AminoAcid M = new("Met", [9, 5, 1, 1, 1], 0.0229);
    public static readonly AminoAcid N = new("Asn", [6, 4, 2, 2, 0], 0.0382);
    public static readonly AminoAcid P = new("Pro", [7, 5, 1, 1, 0], 0.0524);
    public static readonly AminoAcid Q = new("Gln", [8, 5, 2, 2, 0], 0.0394);
    public static readonly AminoAcid R = new("Arg", [12, 6, 4, 1, 0], 0.0590);
    public static readonly AminoAcid S = new("Ser", [5, 3, 1, 2, 0], 0.0733);
    public static readonly AminoAcid T = new("Thr", [7, 4, 1, 2, 0], 0.0557);
    public static readonly AminoAcid V = new("Val", [9, 5, 1, 1, 0], 0.0674);
    public static readonly AminoAcid W = new("Trp", [10, 11, 2, 1, 0], 0.0127);
    public static readonly AminoAcid Y = new("Tyr", [9, 9, 1, 2, 0], 0.0282);
    public static readonly AminoAcid[] All = [A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y];
}
