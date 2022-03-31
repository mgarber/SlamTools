package edu.umms.garberlab.slam;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class TestSplitterWorkings {

	@Test
	void testReadNameMatching() {
		String readName1 = "NS500602:1107:H5V5NBGXL:1:11101:25053:1144 1:N:0:CAACTCTC+GAGAAGAT";
		String readName2 = "NS500602:1107:H5V5NBGXL:1:11101:25053:1144 2:N:0:CAACTCTC+GAGAAGAT";
		
		assertTrue(Hisat3nSplitterWriter.readNamesMatch(readName1, readName2));
		
		readName1 = "NS500602:1107:H5V5NBGXL:1:11101:25053:1144";
		readName2 = "NS500602:1107:H5V5NBGXL:1:11101:25053:1144";
		
		assertTrue(Hisat3nSplitterWriter.readNamesMatch(readName1, readName2));
		
		
	}

}
