// ==========================================================================
// Fast-Matrix-Multiplication library
// Authors: J-G. Dumas, C. Pernet, A. Sedoglavic
// ==========================================================================

/********************************************************************
 * Replacing variable names in Matrix-Vector straight-line programs
 * Reference:
 *   [ J-G. Dumas, C. Pernet, A. Sedoglavic;
 *     Strassen's algorithm is not optimally accurate
 *     ISSAC 2024, Raleigh, NC USA, pp. 254-263.
 *     (https://hal.science/hal-04441653) ]
 ****************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <regex>
#include <cassert>

// Replaces numbered 'input' (from 'length' to '0') to:
//          numbered 'input' appended with 'replace'
std::string& appendchar(std::string& program, const std::string& input,
			const std::string& replace, const size_t length, const char comment) {
    std::clog << comment << " Appendchar  " << input << "0.." << (length-1) << " --> "
	      << input << replace << "0.." << (length-1) << std::endl;
    for(int k(length); --k >= 0;) {
	const std::string source(input+std::to_string(k));
	const std::string target(input+replace+std::to_string(k));
// std::clog << comment << "  replace " << source << " --> " << target << std::endl;
	program = std::regex_replace(program,
				     std::regex(source),
				     target);
    }
    return program;
}


// Replaces 1D-numbered 'input' (from 'length' to '0') to:
//          2D-numbered 'replace', row-major, (s/n) times n
std::string& unvectorize(std::string& program, const std::string& input,
                         const std::string& replace, const size_t n,
                         const size_t s, const char comment) {

    std::clog << comment << " Unvectorize " << input << "0.." << (s-1) << " --> "
	      << replace << "[1.." << (s/n) << ',' << "1.." << n << ']'
	      << std::endl;
    for(int k(s); --k >= 0;) {
	size_t i(k/n), j(k%n);
	const std::string source(input+std::to_string(k));
	std::stringstream target;
	target << replace << '[' << (i+1) << ',' << (j+1) << ']';

//         std::clog << comment << "  replace " << source << " --> " << target.str() << std::endl;
	program = std::regex_replace(program,
				     std::regex(source),
				     target.str());
    }
    return program;
}



// Main Variable Replacer from SLP (input or output) to Matlab Program
// Examples:
//  [+] ./replacer L_3x6x3.slp i o A 3 6 40 1
//      replaces i 0..17 by A(r 0..2,c 0..5)
//      replaces o 0..39 by oA 0..39
//  [+] ./replacer P_3x6x3.slp i o C 3 3 40 0
//      replaces i 0..39 by iC 0..39
//      replaces o 0..8 by oC 0..8
//      produces C = [ oC0..2 ; oC3..5 ; oC6..8 ] ;

int MatlabVariableReplacer(std::ostream& sout, std::istream& SLP,
			   const std::string& inchar, const std::string& ouchar,
			   const std::string rechar, const size_t m, const size_t n,
			   const size_t s, const size_t r, const size_t iotype) {

    const char comment('%');
    std::stringstream buffer; buffer << SLP.rdbuf();
    std::string program(buffer.str());

    std::clog << comment << "  replace := --> =" << std::endl;
    program = std::regex_replace(program, std::regex(":="), " = ");

    if (iotype == 1) {
	const char ichar('r');
	const char jchar('c');
	assert(ichar != inchar[0]);
	assert(jchar != inchar[0]);

	sout << "[m,n] = size(" << rechar << ");" << std::endl;
	sout << "m0 = 0;";
	for(size_t k(1); k<m; ++k) {
	    sout << " m" << k << " = " << k << "*m/" << m << ';';
	}
	sout << " m" << m << " = m;" << std::endl;
	for(size_t k(0); k<m; ++k) {
	    sout << ' ' << ichar << k << " = m" << k << "+1:m" << (k+1) << ';';
	}
	sout << std::endl;

	sout << "n0 = 0;";
	for(size_t k(1); k<n; ++k) {
	    sout << " n" << k << " = " << k << "*n/" << n << ';';
	}
	sout << " n" << n << " = n;" << std::endl;
	for(size_t k(0); k<n; ++k) {
	    sout << ' ' << jchar << k << " = n" << k << "+1:n" << (k+1) << ';';
	}
	sout << std::endl;


	std::clog << comment << " Unvectorize " << inchar << "0.." << (s-1) << " --> "
		  << rechar << "0.." << (m-1) << 'x' << "0.." << (n-1) << std::endl;
	for(int k(s); --k >= 0;) {
	    size_t i(k/n), j(k%n);
	    const std::string source(inchar+std::to_string(k));
	    std::stringstream target;
	    target << rechar << '(' << ichar << i << ',' << jchar << j << ')';

//             std::clog << comment << "  replace " << source << " --> " << target.str()
//                       << std::endl;
	    program = std::regex_replace(program,
					 std::regex(source),
					 target.str());
	}

    }

    if (iotype == 0) {
        appendchar(program, inchar, rechar, r, comment);
    }

    appendchar(program, ouchar, rechar, (iotype==1?r:s), comment);


    sout << program  << std::endl;

    if (iotype == 0) {
	sout << rechar << " = [";
	size_t k(0);
	for(size_t i(0); i<m; ++i) {
	    for(size_t j(0); j<n; ++j, ++k) {
		sout << ' ' << ouchar << rechar << k;
	    }
	    if (i != (m-1)) sout << ' ' << ';';
	}
	sout << " ] ;" << std::endl;
    }

    return 0;
}

// Main Variable Replacer from SLP (input or output) to Maple Program
// Examples:
//  [+] ./replacer L_3x6x3.slp i o A 3 6 40 1
//      replaces i 0..17 by A[1..3, 1..6] (row-major)
//      replaces o 0..39 by oA 0..39
//  [+] ./replacer P_3x6x3.slp i o C 3 3 40 0
//      replaces i 0..39 by iC 0..39
//      replaces o 0..8 by C[1..3,1..3]

int MapleVariableReplacer(std::ostream& sout, std::istream& SLP,
			  const std::string& inchar, const std::string& ouchar,
			  const std::string rechar, const size_t m, const size_t n,
			  const size_t s, const size_t r, const size_t iotype) {

    const char comment('#');
    std::stringstream buffer; buffer << SLP.rdbuf();
    std::string program(buffer.str());

//     if (iotype == 0) {
//         sout << rechar << ":=Matrix(" << m << ',' << n << ");" << std::endl;
//     }

    if (iotype == 1) {
        unvectorize(program, inchar, rechar, n, s, comment);
    }

    appendchar(program, (iotype==1 ? ouchar : inchar), rechar, r, comment);

    if (iotype == 0) {
        unvectorize(program, ouchar, rechar, n, s, comment);
    }

    sout << program  << std::endl;

    return 0;
}

int VariableReplacer(std::ostream& sout, std::istream& SLP,
		     int argc, char **argv, size_t startc) {
    size_t offset(startc);
    bool Matlab(true);
    const std::string opts(argv[offset+1]);
    if (opts == "-M") { ++offset; Matlab = false; }	// Maple
    if (opts == "-L") { ++offset; Matlab = true;  }


    const std::string inchar(argv[++offset]);// input variables common name
    const std::string ouchar(argv[++offset]);// output variables common name
    const std::string rechar(argv[++offset]);// replaced variables common name

    const size_t m(atoi(argv[++offset]));	 // Matrix row dimension
    const size_t n(atoi(argv[++offset]));	 // Matrix column dimension
    const size_t s(m*n);
    const size_t r(atoi(argv[++offset]));	 // Tensor rank
    ++offset;								 // Input (1) or Output (0) matrix
    const size_t iotype(argc>offset?atoi(argv[offset]):0);

    if (Matlab)
	return MatlabVariableReplacer(sout, SLP, inchar, ouchar, rechar, m, n, s, r, iotype);
    else
	return MapleVariableReplacer(sout, SLP, inchar, ouchar, rechar, m, n, s, r, iotype);
}


int main(int argc, char **argv) {
    if ( (argc<7) || (argv[1] == "-h")) {
	    std::clog << "Usage: " << argv[0]
                  << " [stdin|file.slp] [-M|-L] i o A m n r [0|1]\n"
                  << " WARNING, with -L&0: variables m#/n#/r#/c# are reserved (forbidden in SLP)\n"
                  << " -M/-L: maple/matlab output\n"
                  << " i/o/A: input/ouput/additional character\n"
                  << " m/n/r: (mxn, 2D) and (r, 1D) replacements\n"
                  << " 0: o is replaced by A numbered (#mxn), i by iA numbered (#r)\n"
                  << " 1: i is replaced by A numbered (#mxn), o by oA numbered (#r)\n";
	    exit(-1);
    }

    std::ifstream filename(argv[1]);
    if ( filename ) {
        return VariableReplacer(std::cout, filename, argc, argv, 1);
    } else {
        return VariableReplacer(std::cout, std::cin, argc, argv, 0);
    }

}
