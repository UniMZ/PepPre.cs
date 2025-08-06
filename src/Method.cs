namespace UniMZ.PepPre;

using UniMZ.Core;

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
    public static List<(int i, int j, double mz)>[] BuildConstraints(Memory<Peak<double>> peak, Ion[] ion, double e,
        IsotopePatternVector ipv)
    {
        var cs = Enumerable.Range(0, peak.Length + 1).Select(v => new List<(int i, int j, double mz)>()).ToArray();
        for (var i = 0; i < ion.Length; i++)
        {
            var m = ipv.Mass(ion[i].M());
            for (var j = 0; j < m.Length; j++)
            {
                var mz = m[j] / ion[i].Z + ion[i].MZ;
                var idx = peak.ArgQueryNearest(mz);
                if (peak.IsEmpty || Math.Abs(peak.Span[idx].MZ - mz) > e) // empty
                {
                    cs[^1].Add((i, j, mz));
                }
                else
                {
                    cs[idx].Add((i, j, mz));
                }
            }
        }
        return cs;
    }
}
