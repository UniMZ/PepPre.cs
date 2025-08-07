namespace UniMZ.PepPreTest;

using UniMZ.Core;
using UniMZ.PepPre;

[TestClass]
public class MethodTest
{
    [TestMethod]
    public void TestGenerateCandidate()
    {
        var ipv = PepPre.BuildIsotopePatternVector(1000, DefaultIsotopeTable.HCNOS,
            PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All));
        Peak<double>[] peak =
        [
            new(10, 1.0),
            new(100.0, 1.0), new(100.33, 1.0), new(100.5, 1.0), new(100.66, 1.0), new(101, 1.0),
            new(200, 1.0), new(200.5, 1.0)
        ];
        Ion[] ans = [new(100, 2), new(100, 3), new(100.33, 3), new(100.5, 2), new(100.66, 3)];
        var ions = PepPre.GenerateCandidate(peak.AsMemory(), 100, 2, 2, 4, ipv, 100e-6);
        CollectionAssert.AreEqual(ans, ions);
    }

    [TestMethod]
    public void TestBuildConstraints()
    {
        var ipv = PepPre.BuildIsotopePatternVector(1000, DefaultIsotopeTable.HCNOS,
            PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All));
        const double e = 1e-2;
        Peak<double>[] peak =
        [
            new(100, 1.0), new(100.33, 2.0), new(100.330001, 3.0), new(100.66, 4.0)
        ];
        Ion[] ion =
        [
            new(99, 1), new(100, 2), new(100, 3), new(100.33, 3), new(200, 2)
        ];
        var cs = PepPre.BuildConstraints(peak, ion, e, ipv);
        (int i, int j)[][] ans =
        [
            [new ValueTuple<int, int>(0, 1), new ValueTuple<int, int>(1, 0), new ValueTuple<int, int>(2, 0)],
            [new ValueTuple<int, int>(3, 0)],
            [new ValueTuple<int, int>(2, 1)],
            [new ValueTuple<int, int>(2, 2), new ValueTuple<int, int>(3, 1)],
            [
                new ValueTuple<int, int>(0, 0),
                new ValueTuple<int, int>(1, 1), new ValueTuple<int, int>(1, 2), new ValueTuple<int, int>(3, 2),
                new ValueTuple<int, int>(4, 0), new ValueTuple<int, int>(4, 1), new ValueTuple<int, int>(4, 2)
            ]
        ];
        Assert.AreEqual(peak.Length + 1, cs.Length);
        for (var i = 0; i < cs.Length; i++)
        {
            Assert.AreEqual(ans[i].Length, cs[i].Count);
            for (var j = 0; j < cs[i].Count; j++)
            {
                Assert.AreEqual(ans[i][j].i, cs[i][j].i);
                Assert.AreEqual(ans[i][j].j, cs[i][j].j);
                if (i != cs.Length - 1)
                {
                    Assert.IsTrue(Math.Abs(cs[i][j].mz - peak[i].MZ) / peak[i].MZ <= e);
                    Assert.AreEqual(ipv.Abundance(ion[cs[i][j].i].M(), cs[i][j].j), cs[i][j].abu);
                }
            }
        }
    }
}
