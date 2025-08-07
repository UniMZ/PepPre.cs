namespace UniMZ.PepPre;

using UniMZ.Core;
using UniMZ.GLPK;

public static partial class PepPre
{
    public static Ion[] GenerateCandidate(Memory<Peak<double>> peak, double mz, double r, int z_min, int z_max,
        IsotopePatternVector ipv, double error)
    {
        var ions = new Ion[(z_max - z_min + 1) * peak.Length];
        var i = 0;
        peak = peak.Query(mz - r - 2, mz + r + 1);
        foreach (var p in peak.Span)
        {
            for (var z = z_min; z <= z_max; z++)
            {
                var n = (int)(p.MZ * z);
                if (n < ipv.MaxMass() && !peak.QueryError(p.MZ + ipv.Mass(n, 1) / z, error).IsEmpty)
                {
                    ions[i++] = new Ion(p.MZ, z);
                }
            }
        }
        return ions[..i];
    }

    // return an array of lists, where the i-th elements are constraints for the i-th peak, and the last element is for the virtual peak
    public static List<(int i, int j, double mz, double abu)>[] BuildConstraints(Memory<Peak<double>> peak, Ion[] ion,
        double e, IsotopePatternVector ipv)
    {
        var cs = Enumerable.Range(0, peak.Length + 1).Select(v => new List<(int i, int j, double mz, double abu)>())
            .ToArray();
        for (var i = 0; i < ion.Length; i++)
        {
            var m = ipv.Mass(ion[i].M());
            var a = ipv.Abundance(ion[i].M());
            for (var j = 0; j < m.Length; j++)
            {
                var mz = m[j] / ion[i].Z + ion[i].MZ;
                var idx = peak.ArgQueryNearest(mz);
                if (peak.IsEmpty || Math.Abs(peak.Span[idx].MZ - mz) > e) // empty
                {
                    cs[^1].Add((i, j, mz, a[j]));
                }
                else
                {
                    cs[idx].Add((i, j, mz, a[j]));
                }
            }
        }
        return cs;
    }

    public static (double[] Weight, double[] Error) SolveProblem(Memory<Peak<double>> peak, Ion[] ion,
        List<(int i, int j, double mz, double abu)>[] con)
    {
        var nc = ion.Length + con.Length - 1;
        var nr = (con.Length - 1) * 2;
        var prob = GLPK.CreateProb();
        GLPK.AddRows(prob, nr);
        GLPK.AddCols(prob, nc);

        var s = new double[ion.Length];
        foreach (var c in con[^1])
        {
            s[c.i] += c.abu;
        }
        for (var i = 1; i <= ion.Length; i++) GLPK.SetObjCoef(prob, i, s[i - 1]);
        for (var i = ion.Length + 1; i <= nc; i++) GLPK.SetObjCoef(prob, i, 1);
        GLPK.SetObjDir(prob, GLPK.MIN);

        for (var i = 1; i <= nc; i++) GLPK.SetColBnds(prob, i, GLPK.LO, 0, 0);

        for (var i = 0; i < con.Length - 1; i++)
        {
            var len = con[i].Count + 1;
            var idxs = new int[len + 1];
            var vals1 = new double[len + 1];
            var vals2 = new double[len + 1];
            for (var j = 0; j < con[i].Count; j++)
            {
                idxs[j + 1] = con[i][j].i + 1;
                vals1[j + 1] = con[i][j].abu;
                vals2[j + 1] = -con[i][j].abu;
            }
            idxs[^1] = i + 1 + ion.Length;
            vals1[^1] = 1.0;
            vals2[^1] = 1.0;
            GLPK.SetRowBnds(prob, i * 2 + 1, GLPK.LO, peak.Span[i].Inten, 0);
            GLPK.SetRowBnds(prob, i * 2 + 2, GLPK.LO, -peak.Span[i].Inten, 0);
            GLPK.SetMatRow(prob, i * 2 + 1, len, idxs, vals1);
            GLPK.SetMatRow(prob, i * 2 + 2, len, idxs, vals2);
        }

        GLPK.InitSmcp(out var cfg);
        cfg.presolve = GLPK.ON;
        cfg.msg_lev = GLPK.MSG_OFF;
        GLPK.Simplex(prob, ref cfg);

        var ws = Enumerable.Range(1, ion.Length).Select(i => GLPK.GetColPrim(prob, i)).ToArray();
        var es = Enumerable.Range(ion.Length + 1, con.Length - 1).Select(i => GLPK.GetColPrim(prob, i)).ToArray();
        GLPK.DeleteProb(prob);
        return (ws, es);
    }
}
