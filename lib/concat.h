/*
Copyright (C) 2016-present Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
#include <getopt.h>

#include "utility.h"
#include "two_reader.h"

void concat_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Concatenate two or more TWO files\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " concat [options] -i <in.two> -i <in.two> -o <out.two>\n\n"
	"Options:\n"
	"  -i FILE    input TWO file specified 1-or-more times (required)\n"
	"  -I STRING  input file list (required)\n"
	"  -o FILE    output file (- for stdout; default: -)\n" << std::endl;
}

/**<
 * Reads a list of strings and stores them in the provided reference vector
 * of strings. This function is used to grab filenames from lists of strings.
 * @param list  Src vector of file-names pointing to lists of file-names.
 * @param files Dst vector of strings.
 * @return      Returns TRUE upon success or FALSE otherwise.
 */
bool GrabNames(const std::vector<std::string>& list, std::vector<std::string>& files){
	for(int i = 0; i < list.size(); ++i){
		std::ifstream stream(list[i], std::ios::in);
		if(stream.good() == false){
			std::cerr << "faield to open list=" << list[i] << std::endl;
			return false;
		}

		std::string line;
		while(getline(stream,line)){
			//std::cerr << "adding file=" << line << std::endl;
			files.push_back(line);
		}
	}
	return true;
}

int concat(int argc, char** argv){
	if(argc < 3){
		concat_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",   optional_argument, 0, 'i' },
		{"output",  optional_argument, 0, 'o' },
		{"list",    optional_argument, 0, 'I' },

		{0,0,0,0}
	};

	std::vector<std::string> in_list;
	std::vector<std::string> in_file_list;
	std::string out;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:o:?", long_options, &long_index)) != -1){
		hits += 2;
		switch (c){
		case ':':   /* missing option argument */
			fprintf(stderr, "%s: option `-%c' requires an argument\n",
					argv[0], optopt);
			break;

		case '?':
		default:
			fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
					argv[0], optopt);
			break;

		case 'i':
			in_list.push_back(std::string(optarg));
			break;
		case 'o':
			out = std::string(optarg);
			break;
		case 'I':
			in_file_list.push_back(std::string(optarg));
			break;
		}
	}

	if(in_list.size() == 0 && in_file_list.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(in_file_list.size()){
		if(GrabNames(in_file_list, in_list) == false){
			return 1;
		}
	}

	if(in_list.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(in_list.size() == 1){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Only one input file provided..." << std::endl;
		return(1);
	}

	tomahawk::two_reader oreader;

	// open first
	if(oreader.Open(in_list[0]) == false) return 1;

	// Open each file and peek at the headers to check the files are possible to merge.
	for(int i = 1; i < in_list.size(); ++i){
		tomahawk::two_reader rdr;

		// open first
		if(rdr.Open(in_list[i]) == false){
		    std::cerr << tomahawk::utility::timestamp("ERROR") << "Failed to open \"" << in_list[i] << "\"..." << std::endl;
			return 1;
		}

		if(oreader.hdr.samples_.size() != rdr.hdr.samples_.size()){
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Sample have different sample lengths..." << std::endl;
			return 1;
		}

		for(int j = 0; j < oreader.hdr.samples_.size(); ++j){
			if(oreader.hdr.samples_[j] != rdr.hdr.samples_[j]){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Sample have different sample names..." << std::endl;
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Conflict: " << oreader.hdr.samples_[j] << "!=" << rdr.hdr.samples_[i] << " in file " << in_list[i] << std::endl;
				return 1;
			}
		}
	}
	std::cerr << tomahawk::utility::timestamp("LOG") << "All files are compatible. Beginning merging..." << std::endl;

	std::string concat_string = "\n##tomahawk_concatVersion=" + std::string(VERSION) + "\n";
	concat_string += "##tomahawk_concatCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + tomahawk::utility::datetime(); + "\n";
	oreader.hdr.literals_ += concat_string;

	tomahawk::twk_two_writer_t writer;
	writer.mode = 'b';
	writer.oindex.SetChroms(oreader.hdr.GetNumberContigs());

	// Take care of output suffix.
	std::string base_path = tomahawk::twk_writer_t::GetBasePath(out);
	std::string base_name = tomahawk::twk_writer_t::GetBaseName(out);
	std::string extension = tomahawk::twk_writer_t::GetExtension(out);
	if(extension.length() == 3){
		if(strncasecmp(&extension[0], "two", 3) != 0){
			out = (base_path.size() ? base_path + "/" : "") + base_name + ".two";
		}
	} else {
		 out = (base_path.size() ? base_path + "/" : "") + base_name + ".two";
	}

	std::cerr << tomahawk::utility::timestamp("LOG","WRITER") << "Opening " << out << "..." << std::endl;
	if(writer.Open(out) == false){
		std::cerr << "failed to open " << out << std::endl;
		return false;
	}

	writer.WriteHeader(oreader);

	// Loggers.
	uint64_t n_b = 0, n_bc = 0;
	uint64_t nt_b = 0, nt_bc = 0;

	// Begin merge
	// Step0: write out first file
	std::cerr << tomahawk::utility::timestamp("LOG") << "Appending " << in_list[0] << "... ";
	std::cerr.flush();
	uint32_t idx_offset = 0;
	while(oreader.NextBlockRaw()){ // get raw uncompressed data
		tomahawk::IndexEntryOutput rec = oreader.index.ent[idx_offset];
		rec.foff = writer.stream.tellp();
		writer.stream << oreader.it.oblk;
		rec.fend = writer.stream.tellp();
		rec.b_cmp = oreader.it.oblk.nc;
		rec.b_unc = oreader.it.oblk.n;
		n_b  += rec.b_unc;
		n_bc += rec.b_cmp;
		writer.oindex += rec;
		++idx_offset;
	}
	std::cerr << tomahawk::utility::ToPrettyDiskString(n_b) << "/" << tomahawk::utility::ToPrettyDiskString(n_bc) << std::endl;
	nt_b += n_b; nt_bc += n_bc;
	n_b = 0, n_bc = 0;

	for(int i = 1; i < in_list.size(); ++i){
		tomahawk::two_reader rdr;

		// open first
		if(rdr.Open(in_list[i]) == false){
			std::cerr << "failed to open=" << in_list[i] << std::endl;
			return 1;
		}

		std::cerr << tomahawk::utility::timestamp("LOG") << "Appending " << in_list[i] << "... ";
		std::cerr.flush();

		idx_offset = 0;
		while(rdr.NextBlockRaw()){ // get raw uncompressed data
			tomahawk::IndexEntryOutput rec = rdr.index.ent[idx_offset];
			rec.foff = writer.stream.tellp();
			writer.stream << rdr.it.oblk;
			rec.fend = writer.stream.tellp();
			rec.b_cmp = rdr.it.oblk.nc;
			rec.b_unc = rdr.it.oblk.n;
			n_b  += rec.b_unc;
			n_bc += rec.b_cmp;
			writer.oindex += rec;
			++idx_offset;
		}
		std::cerr << tomahawk::utility::ToPrettyDiskString(n_b) << "/" << tomahawk::utility::ToPrettyDiskString(n_bc) << std::endl;
		nt_b += n_b; nt_bc += n_bc;
		n_b = 0, n_bc = 0;
	}

	std::cerr << tomahawk::utility::timestamp("LOG") << "Finished. Added " << in_list.size() << " files..." << std::endl;
	std::cerr << tomahawk::utility::timestamp("LOG") << "Total size: Uncompressed = " << tomahawk::utility::ToPrettyDiskString(nt_b) <<  " and compressed = " << tomahawk::utility::ToPrettyDiskString(nt_bc) << std::endl;

	if(writer.mode == 'b') writer.WriteFinal();
	else writer.WriteBlock();
	writer.close();
	return(0);
}
