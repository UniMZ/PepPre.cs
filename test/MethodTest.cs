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
}
