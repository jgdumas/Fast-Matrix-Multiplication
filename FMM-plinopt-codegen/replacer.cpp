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
std::string& unvectorize(std::string& program, const std::string& input,
                         const std::string& replace, const size_t length) {
    std::clog << "% Unvectorize " << input << "0.." << (length-1) << " --> "
              << input << replace << "0.." << (length-1) << std::endl;
    for(int k(length); --k >= 0;) {
        const std::string source(input+std::to_string(k));
        const std::string target(input+replace+std::to_string(k));
//         std::clog << "%  replace " << source << " --> " << target
//                   << std::endl;
        program = std::regex_replace(program,
                                     std::regex(source),
                                     target);
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

int VariableReplacer(std::ostream& sout, std::istream& SLP,
                     int argc, char **argv, size_t startc) {

    std::stringstream buffer; buffer << SLP.rdbuf();
    std::string program(buffer.str());
    size_t offset(startc);

    const std::string inchar(argv[++offset]);// input variables common name
    const std::string ouchar(argv[++offset]);// output variables common name
    const std::string rechar(argv[++offset]);// replaced variables common name

    const size_t m(atoi(argv[++offset]));	 // Matrix row dimension
    const size_t n(atoi(argv[++offset]));	 // Matrix column dimension
    const size_t s(m*n);
    const size_t r(atoi(argv[++offset]));	 // Tensor rank
    ++offset;								 // Input or Output matrix
    const bool isin(argc>offset?atoi(argv[offset]):true);


    std::clog << "%  replace := --> =" << std::endl;
    program = std::regex_replace(program, std::regex(":="), " = ");

    const char ichar('r');
    const char jchar('c');
    assert(ichar != inchar[0]);
    assert(jchar != inchar[0]);

    if (isin) {
        sout << "[m,n] = size(" << rechar << ");" << std::endl;
        sout << "m0 = 0;";
        for(size_t k(1); k<m; ++k) {
            sout << " m" << k << " = " << k << "*m/" << m << ';';
        }
        sout << " m" << m << " = m;" << std::endl;
        for(size_t k(0); k<m; ++k) {
            sout << ichar << k << " = m" << k << "+1:m" << (k+1) << "; ";
        }
        sout << std::endl;

        sout << "n0 = 0;";
        for(size_t k(1); k<n; ++k) {
            sout << " n" << k << " = " << k << "*n/" << n << ';';
        }
        sout << " n" << n << " = n;" << std::endl;
        for(size_t k(0); k<n; ++k) {
            sout << jchar << k << " = n" << k << "+1:n" << (k+1) << "; ";
        }
        sout << std::endl;

        std::clog << "% Unvectorize " << inchar << "0.." << (s-1) << " --> "
                  << rechar << "0.." << (m-1) << 'x' << "0.." << (n-1) << std::endl;
        for(int k(s); --k >= 0;) {
            size_t i(k/n), j(k%n);
            const std::string source(inchar+std::to_string(k));
            std::stringstream target;
            target << rechar << '(' << ichar << i << ',' << jchar << j << ')';

//             std::clog << "%  replace " << source << " --> " << target.str()
//                       << std::endl;
            program = std::regex_replace(program,
                                         std::regex(source),
                                         target.str());
        }

    } else {
        unvectorize(program, inchar, rechar, r);
    }

    unvectorize(program, ouchar, rechar, (isin?r:s) );


    sout << program  << std::endl;

    if (! isin) {
        sout << "C = [";
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



int main(int argc, char **argv) {
    if ( (argc<7) || (argv[1] == "-h")) {
            std::clog << "Usage: " << argv[0]
                      << " [stdin|file.slp] i o A m n r [0|1]"
                      << std::endl;
            exit(-1);
    }

    std::ifstream filename(argv[1]);
    if ( filename ) {
        return VariableReplacer(std::cout, filename, argc, argv, 1);
    } else {
        return VariableReplacer(std::cout, std::cin, argc, argv, 0);
    }

}
