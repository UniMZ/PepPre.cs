namespace UniMZ.PepPre;

using System.Numerics;
using UniMZ.Core;

public readonly record struct IsotopePatternVector(double[][] Abundance, double[][] Mass, int[] MaxIndex);

public static partial class PepPre
{
    public static (TW[] Weight, TV[] Value) Conv<TW, TV>(TW[] wa, TW[] wb, TV[] va, TV[] vb)
        where TW : INumber<TW>, IMultiplyOperators<TV, TW, TV>
        where TV : INumber<TV>, IMultiplyOperators<TV, TW, TV>, IDivisionOperators<TV, TW, TV>
    {
        var len = wa.Length + wb.Length - 1;
        var w = new TW[len];
        var v = new TV[len];
        for (var i = 0; i < len; i++)
        {
            var l = Math.Max(0, i - wb.Length + 1);
            var n = Math.Min(i + 1, wa.Length) - l;
            var ws = Enumerable.Range(l, n).Select(x => wa[x] * wb[i - x]).ToArray();
            var vs = Enumerable.Range(l, n).Select(x => va[x] + vb[i - x]).ToArray();
            w[i] = ws.Sum();
            v[i] = w[i] != default ? vs.Mul(ws).Sum() / w[i] : default!;
        }
        return (w, v);
    }

    public static double[] AverageFormulaPerAminoAcid(AminoAcid[] tab)
    {
        var n = tab[0].Formula.Length;
        var s = tab.Select(v => v.Abundance).Sum();
        return Enumerable.Range(0, n).Select(i => tab.Select(v => v.Formula[i] * v.Abundance).Sum() / s).ToArray();
    }

    public static double[] AverageFormulaPerDalton(Isotope[] tab, double[] formula)
    {
        return formula.Div(tab.Select((v, i) => v.Mass[0] * formula[i]).Sum()).ToArray();
    }

    public static double AverageMassPerAminoAcid(AminoAcid[] aa, Isotope[] iso)
    {
        var x = AverageFormulaPerAminoAcid(aa);
        return iso.Select((v, i) => v.Mass[0] * x[i]).Sum();
    }

    // make it incremental
    public static void SmoothMaxIndex(this IsotopePatternVector ipv)
    {
        for (var i = 1; i < ipv.MaxIndex.Length; ++i) ipv.MaxIndex[i] = Math.Max(ipv.MaxIndex[i - 1], ipv.MaxIndex[i]);
    }

    public static IsotopePatternVector BuildIsotopePatternVector(int max_mass, Isotope[] tab, double[] formula,
        double trunc = 0.99, bool smooth = true)
    {
        var abundance = new double[max_mass + 1][];
        var mass = new double[max_mass + 1][];
        var idx = new int[max_mass + 1];
        abundance[0] = [1.0];
        mass[0] = [0.0];
        idx[0] = 0;
        var n_mean = AverageFormulaPerDalton(tab, formula);
        var n_max = n_mean.Select(v => (int)Math.Round(v * max_mass)).ToArray();
        var ws = n_max.Select(v => new double[Math.Max(2, v + 1)][]).ToArray(); // weight
        var vs = n_max.Select(v => new double[Math.Max(2, v + 1)][]).ToArray(); // value
        for (var i = 0; i < tab.Length; i++)
        {
            ws[i][0] = [1.0];
            vs[i][0] = [0.0];
            var sum = tab[i].Abundance.Sum();
            var mono = tab[i].Mass[0];
            ws[i][1] = tab[i].Abundance.Select(v => v / sum).ToArray();
            vs[i][1] = tab[i].Mass.Select(v => v - mono).ToArray();
            for (var j = 2; j < ws[i].Length; j++)
                (ws[i][j], vs[i][j]) = Conv(ws[i][1], ws[i][j - 1], vs[i][1], vs[i][j - 1]);
        }
        Parallel.For(1, max_mass + 1, m =>
        {
            var ns = n_mean.Select(v => Math.Round(v * m)).ToArray();
            var w = ws[0][(int)ns[0]];
            var v = vs[0][(int)ns[0]];
            for (var i = 1; i < ns.Length; i++)
            {
                if (ns[i] < 0) continue;
                (w, v) = Conv(w, ws[i][(int)ns[i]], v, vs[i][(int)ns[i]]);
                w = w.Div(w.Sum()).ToArray();
            }
            var s = 0.0;
            var c = 0;
            while (s < trunc)
            {
                s += w[c];
                c += 1;
            }
            abundance[m] = w[..c];
            mass[m] = v[..c];
            idx[m] = abundance[m].Select((value, index) => (value, index)).Max().index;
        });
        var ipv = new IsotopePatternVector(abundance, mass, idx);
        if (smooth) SmoothMaxIndex(ipv);
        return ipv;
    }

    public static IsotopePatternVector ReadIsotopePatternVector(this string path)
    {
        using var fs = File.OpenRead(path);
        using var reader = new BinaryReader(fs);
        var abundance = reader.ReadJaggedArray<double>();
        var mass = reader.ReadJaggedArray<double>();
        var idx = reader.ReadArray<ulong>().Select(v => (int)v - 1).ToArray();
        return new IsotopePatternVector(abundance, mass, idx);
    }

    public static void Write(this IsotopePatternVector ipv, string path)
    {
        var fs = File.Open(path + "~", FileMode.Create);
        var writer = new BinaryWriter(fs);
        writer.WriteJaggedArray(ipv.Abundance);
        writer.WriteJaggedArray(ipv.Mass);
        writer.WriteArray(ipv.MaxIndex.Select(v => (ulong)v + 1).ToArray());
        writer.Close();
        fs.Close();
        File.Delete(path);
        File.Move(path + "~", path);
    }

    public static int MaxMass(this IsotopePatternVector ipv) { return ipv.Mass.Length; }
    public static double[] Mass(this IsotopePatternVector ipv, int m) { return ipv.Mass[m]; }
    public static double[] Mass(this IsotopePatternVector ipv, float m) { return ipv.Mass[(int)m]; }
    public static double[] Mass(this IsotopePatternVector ipv, double m) { return ipv.Mass[(int)m]; }
    public static double[] Abundance(this IsotopePatternVector ipv, int m) { return ipv.Abundance[m]; }
    public static double[] Abundance(this IsotopePatternVector ipv, float m) { return ipv.Abundance[(int)m]; }
    public static double[] Abundance(this IsotopePatternVector ipv, double m) { return ipv.Abundance[(int)m]; }
    public static double Mass(this IsotopePatternVector ipv, int m, int i) { return ipv.Mass[m][i]; }
    public static double Mass(this IsotopePatternVector ipv, float m, int i) { return ipv.Mass[(int)m][i]; }
    public static double Mass(this IsotopePatternVector ipv, double m, int i) { return ipv.Mass[(int)m][i]; }
    public static double Abundance(this IsotopePatternVector ipv, int m, int i) { return ipv.Abundance[m][i]; }
    public static double Abundance(this IsotopePatternVector ipv, float m, int i) { return ipv.Abundance[(int)m][i]; }
    public static double Abundance(this IsotopePatternVector ipv, double m, int i) { return ipv.Abundance[(int)m][i]; }
}
