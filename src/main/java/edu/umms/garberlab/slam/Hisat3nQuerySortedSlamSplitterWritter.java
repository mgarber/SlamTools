package edu.umms.garberlab.slam;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class Hisat3nQuerySortedSlamSplitterWritter extends Hisat3nSplitterWriter {
	private String unconvertedFilePath;
	private String convertedFilePath;
	
	private SAMFileWriter unconvertedWriter;
	private SAMFileWriter convertedWriter;
	
	private List<SAMRecord> currentPairList;
	
	//TODO: either makes this configurable or deal with reads where the name of the pairs is not the same.
	// This can be due to the pair information tacked on to the read name such us in these:
	//NS500602:968:H7FMMBGXC:1:11101:1041:15222 1:N:0:AAGTCCAA        77      *       0       0       *       *       0       0       GCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTCCGGCGCGCCCTGAGACGCGCAGA        AAAAAEEEEEEEEEEAEEEEEEEEE<AEEEE
	//NS500602:968:H7FMMBGXC:1:11101:1041:15222 2:N:0:AAGTCCAA        141     *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	//
	
	public Hisat3nQuerySortedSlamSplitterWritter(File outDirFile, String outPrefix, SamReader reader, boolean asBAM) {
		this.unconvertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_unconverted.bam";
		this.convertedFilePath = outDirFile.getAbsolutePath() + "/" + outPrefix + "_converted.bam";
		
		currentPairList = new ArrayList<SAMRecord>();
		
		SAMProgramRecord pg = new SAMProgramRecord(SlamSplitter.PROGRAM_NAME);
		SAMFileHeader header = reader.getFileHeader();
		header.addProgramRecord(pg);
		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		if(asBAM) {
			unconvertedWriter = factory.makeBAMWriter(header, false, new File(unconvertedFilePath)); 
			convertedWriter = factory.makeBAMWriter(header, false, new File(convertedFilePath));
                    
		} else {
			unconvertedWriter = factory.makeSAMWriter(header, false, new File(unconvertedFilePath));
			convertedWriter = factory.makeSAMWriter(header, false, new File(convertedFilePath));
		}
		
		
	}

	@Override
	public void close() {
		if(!currentPairList.isEmpty()) {
			processPairList();
		}
		unconvertedWriter.close();
		convertedWriter.close();
 
	}

	public void write(SAMRecord samRecord) {
		//TODO: either makes this configurable or deal with reads where the name of the pairs is not the same.
		// This can be due to the pair information tacked on to the read name such us in these:
		//NS500602:968:H7FMMBGXC:1:11101:1041:15222 1:N:0:AAGTCCAA        77      *       0       0       *       *       0       0       GCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTCCGGCGCGCCCTGAGACGCGCAGA        AAAAAEEEEEEEEEEAEEEEEEEEE<AEEEE
		//NS500602:968:H7FMMBGXC:1:11101:1041:15222 2:N:0:AAGTCCAA        141     *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
		
		if(!isWriteUnmapped() && samRecord.getMateUnmappedFlag() && samRecord.getReadUnmappedFlag()) {
			//currentPairList.clear();
			return; // Just ignore them
		}
		
		if(currentPairList.isEmpty() || readNamesMatch(currentPairList.get(0), samRecord)) {
			currentPairList.add(samRecord);
		} else {
			processPairList();
			currentPairList.clear();
			currentPairList.add(samRecord);
		}

	}

	

	public void write(Hisat3nAlignedFragment hisatAlignmentFragment) {
		if(!isPairedEnd() ||  hisatAlignmentFragment.isUnpaired()) {
			write(hisatAlignmentFragment.getPair1(), hisatAlignmentFragment.getConvertedBases()>0);
		} else {
			int convertedBases = hisatAlignmentFragment.getConvertedBases();
			Hisat3nSAMRecord upstreamRecord = hisatAlignmentFragment.getUpstreamRecord();
			write(upstreamRecord, convertedBases > 0);
			Hisat3nSAMRecord downstreamRecord = hisatAlignmentFragment.getDownstreamRecord();
			write(downstreamRecord, convertedBases > 0);
		}

	}
	
	private void write(Hisat3nSAMRecord record, boolean hasConvertedBases) {
		if(hasConvertedBases) {
			if(revertConversion()) {
				record.revertConvertedBases();
			}
			convertedWriter.addAlignment(record.getSAMRecord());
		} else {
			unconvertedWriter.addAlignment(record.getSAMRecord());
		}
	}
	
	
	private void processPairList() {
		if(currentPairList.size() % 2 != 0) {
			System.err.println("ERROR: Found a list that is supposed to contain read pairs with odd size: " + currentPairList.size() + " reads: ");
			for (SAMRecord r : currentPairList) {
				System.err.println(r.toString());
			}
			System.exit(1);
		}
		
		int i = 0;
		while (i < currentPairList.size()) {
			Hisat3nAlignedFragment frag = new Hisat3nAlignedFragment(currentPairList.get(i), currentPairList.get(i+1));
			i = i+2;
			write(frag);
		}
		
	}


}
