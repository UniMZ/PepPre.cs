namespace UniMZ.PepPreTest;

using UniMZ.PepPre;

[TestClass]
public class IsotopePatternVectorTest
{
    [TestMethod]
    public void TestConv()
    {
        double[] wa = [1.1f, 2.2f, 3.3f, 4.4f];
        double[] wb = [2.0f, 3.0f];
        double[] va = [3, 5, 7, 11];
        double[] vb = [2, 3];
        double[] w = [2.2, 7.7, 13.2, 18.7, 13.2];
        double[] v = [5.0, 6.571428571428572, 8.5, 11.411764705882351, 14.0];
        var (wc, vc) = PepPre.Conv(wa, wb, va, vb);
        Assert.IsTrue(w.Zip(wc, (a, b) => Math.Abs(a - b)).Sum() < 1e-6);
        Assert.IsTrue(v.Zip(vc, (a, b) => Math.Abs(a - b)).Sum() < 1e-6);
    }

    [TestMethod]
    public void TestAverageFormulaPerAminoAcid()
    {
        double[] ans = [7.7526, 4.9114, 1.3643, 1.4718, 0.0375];
        var x = PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All);
        for (var i = 0; i < ans.Length; ++i)
            Assert.IsTrue(Math.Abs(x[i] - ans[i]) < 1e-4);
    }

    [TestMethod]
    public void TestAverageFormulaPerDalton()
    {
        double[] ans = [0.07010020, 0.04440957, 0.01233584, 0.01330777, 0.00033905];
        var x = PepPre.AverageFormulaPerDalton(DefaultIsotopeTable.HCNOS,
            PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All));
        for (var i = 0; i < ans.Length; ++i)
            Assert.IsTrue(Math.Abs(x[i] - ans[i]) < 1e-6);
    }

    [TestMethod]
    public void TestAverageMassPerAminoAcid()
    {
        const double ans = 110.59347293;
        var m = PepPre.AverageMassPerAminoAcid(DefaultAminoAcidTable.All, DefaultIsotopeTable.HCNOS);
        Assert.IsTrue(Math.Abs(ans - m) < 1e-4);
    }

    private static void TestIsotopePatternVector(IsotopePatternVector ipv)
    {
        Assert.IsTrue(Math.Abs(ipv.Abundance[1000][0] - 0.57) < 0.01);
        Assert.IsTrue(Math.Abs(ipv.Abundance[1000][1] - 0.30) < 0.01);
        Assert.IsTrue(Math.Abs(ipv.Abundance[1000][2] - 0.09) < 0.01);
        Assert.IsTrue(Math.Abs(ipv.Abundance[1000][3] - 0.02) < 0.01);
        Assert.IsTrue(Math.Abs(ipv.Mass[1000][0] - 0) < 0.0001);
        Assert.IsTrue(Math.Abs(ipv.Mass[1000][1] - 1.0029) < 0.0001);
        Assert.IsTrue(Math.Abs(ipv.Mass[1000][2] - 2.0055) < 0.0001);
        Assert.IsTrue(Math.Abs(ipv.Mass[1000][3] - 3.0080) < 0.0001);
        Assert.IsTrue(ipv.MaxIndex[1000] == 0);
        Assert.IsTrue(ipv.MaxIndex[2000] == 1);
    }

    [TestMethod]
    public void TestBuildIsotopePatternVector()
    {
        var ipv = PepPre.BuildIsotopePatternVector(2000, DefaultIsotopeTable.HCNOS,
            PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All));
        TestIsotopePatternVector(ipv);
    }

    [TestMethod]
    public void TestReadWriteIsotopePatternVector()
    {
        var ipv = PepPre.BuildIsotopePatternVector(2000, DefaultIsotopeTable.HCNOS,
            PepPre.AverageFormulaPerAminoAcid(DefaultAminoAcidTable.All));
        const string path = "tmp.ipv";
        ipv.Write(path);
        var ipv_new = path.ReadIsotopePatternVector();
        TestIsotopePatternVector(ipv_new);
    }
}
