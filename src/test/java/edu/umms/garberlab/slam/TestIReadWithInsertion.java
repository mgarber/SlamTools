package edu.umms.garberlab.slam;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.URL;

public class TestIReadWithInsertion {
	
	
	@Test
	public  void TestInsertionHandling() throws IOException {
		URL inputBamResource = this.getClass().getResource("/insertion.bug.reads.bam");
		assertNotNull(inputBamResource);

		SamInputResource sir = SamInputResource.of(inputBamResource);		
		
		try (SamReader reader1 = SamReaderFactory.makeDefault().open(sir)) {
			try (SamReader reader2 = SamReaderFactory.makeDefault().open(sir)){
				SAMRecordIterator it = reader1.iterator();
				assertTrue(it.hasNext());
				Hisat3nSAMRecord hisatRecord = new Hisat3nSAMRecord(it.next());
				Hisat3nAlignedFragment frag = new Hisat3nAlignedFragment(hisatRecord, reader2);
				assertNotNull(frag);
			}
		}
		
	}

}
