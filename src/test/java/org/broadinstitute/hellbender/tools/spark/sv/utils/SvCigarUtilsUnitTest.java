package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;


public class SvCigarUtilsUnitTest {

    @Test(groups = "sv")
    public void testClippingArithmetic() {
        Cigar cigar = TextCigarCodec.decode("100M51S");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51S100M");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("100M51H");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 0);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 51);

        cigar = TextCigarCodec.decode("51H100M");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 51);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 0);

        cigar = TextCigarCodec.decode("12H12S101M13S13H");
        Assert.assertEquals(SvCigarUtils.getTotalHardClipping(cigar), 25);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(true, cigar), 24);
        Assert.assertEquals(SvCigarUtils.getNumClippedBases(false, cigar), 26);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_emptyCigarElementList(){
        @SuppressWarnings("unchecked")
        final List<CigarElement> emptyList = Collections.EMPTY_LIST;
        SvCigarUtils.validateCigar(emptyList);
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_deletionNeighboringClipping(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10D10M").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S").getCigarElements());
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10M10D10S10H").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_only1NonAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10S").getCigarElements());
    }

    @Test(expectedExceptions=IllegalArgumentException.class, groups = "sv")
    public void testCigarChecker_noAlignment(){
        SvCigarUtils.validateCigar(TextCigarCodec.decode("10H10S10I10S10H").getCigarElements());
    }

    @Test(groups = "sv")
    public void testGetNumClippingBases_hardAndSoftSeparately() {
        List<CigarElement> cigarElements = TextCigarCodec.decode("10H20S30M40D50M60S70H").getCigarElements();
        Assert.assertEquals(SvCigarUtils.getNumHardClippingBases(true, cigarElements), 10);
        Assert.assertEquals(SvCigarUtils.getNumHardClippingBases(false, cigarElements), 70);
        Assert.assertEquals(SvCigarUtils.getNumSoftClippingBases(true, cigarElements), 20);
        Assert.assertEquals(SvCigarUtils.getNumSoftClippingBases(false, cigarElements), 60);
    }

    @Test(groups = "sv")
    public void testGetIndexOfFirstNonClippingBase(){
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), true), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("151M").getCigarElements(), false), 0);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10S10D10M").getCigarElements(), true), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10H10S10D10M").getCigarElements(), true), 2);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S").getCigarElements(), false), 1);
        Assert.assertEquals(SvCigarUtils.findIndexOfFirstNonClippingOperation(TextCigarCodec.decode("10M10D10S10H").getCigarElements(), false), 1);
    }

    @DataProvider(name = "CigarTestData")
    private Object[][] createCigarTestData() {
        final Cigar[] cigars = Stream.of("5M5H", "5H5M", "14S6M", "14M6S",  "12H12S20M51I30M13S13H", "10H20S30M40D50M60S70H").map(TextCigarCodec::decode).toArray(Cigar[]::new);
        final int[] lengths = {10, 10, 20, 20, 151, 240};
        final Object[][] data = new Object[cigars.length][];
        for (int i = 0; i < cigars.length; ++i) {
            data[i] = new Object[]{cigars[i], lengths[i]};
        }
        return data;
    }

    @Test(dataProvider = "CigarTestData", groups = "sv")
    public void testGetUnclippedReadLengthFromCigar(final Cigar cigar, final int expectedReadLength) {
        Assert.assertEquals(SvCigarUtils.getUnclippedReadLength(cigar), expectedReadLength);
    }

    @DataProvider(name = "refWalkDistanceTestDataException")
    private Object[][] createRefWalkDistanceTestDataException() {
        final Object[][] data = new Object[7][];
        data[0] = new Object[]{TextCigarCodec.decode("50M10N101M"), 41, 10, 0};
        data[1] = new Object[]{TextCigarCodec.decode("50M10P101M"), 41, 10, 0};
        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");
        data[2] = new Object[]{cigar, -1, 10, 0};
        data[3] = new Object[]{cigar, 0, 10, 0};
        data[4] = new Object[]{cigar, 41, -1, 0};
        data[5] = new Object[]{cigar, 41, 0, 0};
        data[6] = new Object[]{cigar, 1, 201, 0};
        return data;
    }

    @Test(dataProvider = "refWalkDistanceTestDataException", groups = "sv", expectedExceptions = IllegalArgumentException.class)
    public void testRefWalkDistanceException(final Cigar cigar, final int startInclusive, final int distance,
                                             final int expectedRefDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRef(cigar, startInclusive, distance), expectedRefDist);
    }

    @DataProvider(name = "refWalkDistanceTestData")
    private Object[][] createRefWalkDistanceTestData() {
        final Object[][] data = new Object[11][];
        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");
        data[0] = new Object[]{cigar, 1, 40, 0};
        data[1] = new Object[]{cigar, 1, 45, 5};
        data[2] = new Object[]{cigar, 41, 10, 10};
        data[3] = new Object[]{cigar, 41, 30, 10};
        data[4] = new Object[]{cigar, 41, 25, 10};
        data[5] = new Object[]{cigar, 41, 35, 15};
        data[6] = new Object[]{cigar, 41, 56, 66};
        data[7] = new Object[]{cigar, 41, 110, 115};
        data[8] = new Object[]{cigar, 1, 200, 115};
        data[9] = new Object[]{cigar, 61, 10, 0};
        data[10] = new Object[]{cigar, 61, 15, 5};
        return data;
    }

    @Test(dataProvider = "refWalkDistanceTestData", groups = "sv")
    public void testRefWalkDistance(final Cigar cigar, final int startInclusive, final int distance, final int expectedRefDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRef(cigar, startInclusive, distance), expectedRefDist);
    }

    @DataProvider(name = "readWalkDistanceTestDataException")
    private Object[][] createReadWalkDistanceTestDataException() {
        final Object[][] data = new Object[6][];
        data[0] = new Object[]{TextCigarCodec.decode("50M10N101M"), 41, 10, false, 0};
        data[1] = new Object[]{TextCigarCodec.decode("50M10P101M"), 41, 10, false, 0};
        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");
        data[2] = new Object[]{cigar, -1, 10, false, 0};
        data[3] = new Object[]{cigar, 0, 10, false, 0};
        data[4] = new Object[]{cigar, 41, -1, false, 0};
        data[5] = new Object[]{cigar, 41, 0, false, 0};
        return data;
    }

    @Test(dataProvider = "readWalkDistanceTestDataException", groups = "sv", expectedExceptions = IllegalArgumentException.class)
    public void testReadWalkDistanceTestDataException(final Cigar cigar, final int startInclusive, final int refWalkDist,
                                                      final boolean walkBackwards, final int expectedReadWalkDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRead(cigar, startInclusive, refWalkDist, walkBackwards), expectedReadWalkDist);
    }

    @DataProvider(name = "readWalkDistanceTestData")
    private Object[][] createReadWalkDistanceTestData() {
        final Object[][] data = new Object[4][];
        final Cigar cigar = TextCigarCodec.decode("35H40S10M20I25M30D50M55S60H");
        data[0] = new Object[]{cigar, 1, 10, false, 50};
        data[1] = new Object[]{cigar, 1, 5, false, 45};
        data[2] = new Object[]{cigar, 1, 16, false, 76};
        data[3] = new Object[]{cigar, 1, 66, false, 96};

        return data;
    }

    @Test(dataProvider = "readWalkDistanceTestData", groups = "sv")
    public void testReadWalkDistanceTestData(final Cigar cigar, final int startInclusive, final int refWalkDist,
                                             final boolean walkBackwards, final int expectedReadWalkDist) {
        Assert.assertEquals(SvCigarUtils.computeAssociatedDistOnRead(cigar, startInclusive, refWalkDist, walkBackwards), expectedReadWalkDist);
    }
}
