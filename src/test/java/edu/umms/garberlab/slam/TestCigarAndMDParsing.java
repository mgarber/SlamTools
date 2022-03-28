package edu.umms.garberlab.slam;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

class TestCigarAndMDParsing {
	
	private static final Pattern CigarBlockPattern = Pattern.compile("[0-9]+[S|M|N|I|D]+");

	@Test
	void testLessThanEqualCase() {
		String read = "GTTCCCTTAGCAAGGCTGAAAATTTCAGTCTCTGGTATTTGGAATTTAGGCTGCAGTCCTTGTTTTTGGATGGATCACTG";
		String cigarString = "22M1I27M1I29M";
		String mdTag       = "20G25G1A29";
		String reference   = "GTTCCCTTAGCAAGGCTGAAGA-TTCAGTCTCTGGTATTTGGAATTTGGA-TGCAGTCCTTGTTTTTGGATGGATCACTG";
		
		Cigar cigar = buildCigarFromString(cigarString);
		
		assertEquals(cigarString, cigar.toString(), "Cigar object and cigar string are not the same");
		
		Map<Integer, Character> mismatchMap = Hisat3nSAMRecord.makeSequenceToReferenceMap(cigar, mdTag);
		
		assertTrue(mismatchMap.containsKey(20),"G-->A mismatch at position 20 was not even there");
		assertEquals(mismatchMap.get(20),'G', "Mismatch at position 20 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(22),"Insertion at position 22 was not recorded ");
		assertEquals(mismatchMap.get(22),'-', "Position 22 was recorded as a mismatch but not as an insertion");
		
		assertTrue(mismatchMap.containsKey(47),"G-->A mismatch at position 47 was not even there");
		assertEquals(mismatchMap.get(47),'G', "Mismatch at position 47 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(49),"A-->G mismatch at position 49 was not even there");
		assertEquals(mismatchMap.get(49),'A', "Mismatch at position 49 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(50),"Insertion at position 50 was not recorded ");
		assertEquals(mismatchMap.get(50),'-', "Position 50 was recorded as a mismatch but not as an insertion");
		
		System.out.println(mismatchMap.toString());
		
	}
	
	@Test
	void testMidRefInsertionInMD() {

		String cigarString = "26M1D48M";
		String mdTag       = "16T9^A48";

		
		Cigar cigar = buildCigarFromString(cigarString);
		
		assertEquals(cigarString, cigar.toString(), "Cigar object and cigar string are not the same");
		
		Map<Integer, Character> mismatchMap = Hisat3nSAMRecord.makeSequenceToReferenceMap(cigar, mdTag);
		
		System.out.println(mismatchMap.toString());
		
	}
	
	@Test
	void testTheLessThanPeskyBug() {

		String cigarString = "40M1I39M";
		String mdTag       = "34T5T0C1N34T0";

		
		Cigar cigar = buildCigarFromString(cigarString);
		
		assertEquals(cigarString, cigar.toString(), "Cigar object and cigar string are not the same");
		
		Map<Integer, Character> mismatchMap = Hisat3nSAMRecord.makeSequenceToReferenceMap(cigar, mdTag);
		
		System.out.println(mismatchMap.toString());
		
		assertTrue(mismatchMap.containsKey(34),"T-->C mismatch at position 34 was not even there");
		assertEquals(mismatchMap.get(34),'T', "Mismatch at position 34 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(40),"Insertion at position 40 was not recorded ");
		assertEquals(mismatchMap.get(40),'-', "Position 40 was recorded as a mismatch but not as an insertion");
		
		assertTrue(mismatchMap.containsKey(41),"T-->C mismatch at position 41 was not even there");
		assertEquals(mismatchMap.get(41),'T', "Mismatch at position 41 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(42),"C-->T mismatch at position 42 was not even there");
		assertEquals(mismatchMap.get(42),'C', "Mismatch at position 42 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(44),"N->V mismatch at position 42 was not even there");
		assertEquals(mismatchMap.get(44),'N', "Mismatch at position 41 was not what was expected");
		
		assertTrue(mismatchMap.containsKey(79),"T-->C mismatch at position 79 was not even there");
		assertEquals(mismatchMap.get(79),'T', "Mismatch at position 79 was not what was expected");
		
		
		
	}
	
	
	
	
	private static Cigar buildCigarFromString(String cigarStr) {
		Matcher matcher = CigarBlockPattern.matcher(cigarStr);
		List<CigarElement> elements = new ArrayList<CigarElement>();
		while (matcher.find()) {
			String elementStr = matcher.group();
			int elementLength = Integer.parseInt(elementStr.substring(0, elementStr.length() - 1));
			char opChr        = elementStr.charAt(elementStr.length() - 1);
			CigarOperator co = null;
			
			switch (opChr) {
			case 'S':
				co = CigarOperator.S;
				break;
			case 'M':
				co = CigarOperator.M;
				break;
			case 'N':
				co = CigarOperator.N;
				break;
			case 'I':
				co=CigarOperator.I;
				break;	
			case 'D':
				co=CigarOperator.D;
				break;	
			}
				
				
			
			elements.add(new CigarElement(elementLength, co));
			
		}
		
		return new Cigar(elements);
	}

}
