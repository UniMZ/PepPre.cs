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
}
