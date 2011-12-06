#! /usr/bin/perl
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
$INCLUDE_strg = "INCLUDE =";
########## CHECK FOR NETCDF #############################################
$netcdf_exists = `which ncdump |wc -l `;
if ( $netcdf_exists>0){
  if ( -e nc_routines_true.f90 ) {
  system("mv nc_routines_true.f90 nc_routines.f90");
  }
  print "netcdf exists, hence we will create Makefile for allowing its usage.\n";
  $libs  = "LIBS = -lm -lnetcdf \n";
}
else{
  if ( -e nc_routines_ghost.f90 ) {
  system("mv nc_routines_ghost.f90 nc_routines.f90");
  }
  print "netcdf doesnt exist, hence we will create Makefile disallowing its usage.\n";
  $libs = "LIBS = -lm \n" ;
}
########## CHECK FOR NETCDF #############################################

open(MAKEFILE, "> Makefile");
print MAKEFILE "PROG =xsem\n\n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.f90 *.f *.F *.c>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
print MAKEFILE $libs;
print MAKEFILE "CC = gcc\n";
print MAKEFILE "CFLAGS = -O3 -DF_UNDERSCORE\n";

############ CHOOSE BETWEEN DIFFERENT FORTRAN COMPILERS ###########################
if ($ARGV[0] eq 'ifort'){
    if ($ARGV[1] eq 'debug'){
	$F90_strg = 'mpif90  -vec-report:0 -g -O2 -shared-intel  -mcmodel=medium -ftz -check all -debug  -check -traceback';
	$FC_strg = 'ifort  -vec-report:0 -g -O2 -shared-intel  -mcmodel=medium -ftz -check all -debug  -check -traceback';
    } else {
	$F90_strg = 'mpif90  -vec-report:0 -g -O4 -xHOST -shared-intel'; 
	$FC_strg = 'ifort  -vec-report:0 -g -O4 -xHOST -shared-intel'; 
    }
}
elsif ($ARGV[0] eq 'gfortran'){
    if ($ARGV[1] eq 'debug'){
	$F90_strg = 'mpif90 -Warray-temporaries -fcheck-array-temporaries -fbounds-check -frange-check -pedantic';
	$FC_strg =  'gfortran -Warray-temporaries -fcheck-array-temporaries -fbounds-check -frange-check -pedantic';
    } else {
	$F90_strg = 'mpif90 -O3';
	$FC_strg = 'gfortran  -O3';	
    }
} 
elsif ($ARGV[0] eq '-h'){
	        print "-----------Flags to be used---------- \n";
		        print "Argument 1: Compiler options: gfortran (default), ifort\n";
			        print "Argument 2: debug\n";
				        print "Not specifying debug will create Makefile for optimized compilation \n";
					        exit;
					}
else {
	print "Default compiler is gfortran\n";
	print "If you want another compiler type ./makemake.pl <compiler_name> \n";
	$F90_strg = 'mpif90';
	$FC_strg = 'gfortran';
}
###################################################################################

$F90_full="F90 = $F90_strg \n";
$FC_full="FC = $FC_strg \n";
$INCLUDE_full="INCLUDE = -I $ENV{'NETCDFHOME'}/include \n";

print MAKEFILE $F90_full;
print MAKEFILE $FC_full;
print MAKEFILE $INCLUDE_full; 
print "\n:::::: F90 compiler & flags ::::::\n $F90_strg \n";
print "\n:::::: FC compiler & flags  ::::::\n $FC_strg \n";

print MAKEFILE "\n\n";
print MAKEFILE "# cancel m2c implicit rule \n";
print MAKEFILE "%.o : %.mod \n ";
print MAKEFILE "\n\n";

#
# make
#
print MAKEFILE "all: \$(PROG)\n\n";
print MAKEFILE "\$(PROG): \$(OBJS)\n";
print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.M *.mod *.d *.il core \n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$(INCLUDE) \$<\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs,$1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f90$/.o/;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F90 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F90 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.f90/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }

   print "\nCheck Makefile to make sure you're happy with it.\n";

if ($netcdf_exists>0) {
print "........you may want to double-check the netcdf library path in the makefile...\n";
print "........you may need to include further netcdf libraries (libnetcdf_c++,libnetcdff,lcurl)\n";
print "........you may wish to find the appropriate INCLUDE (locate NETCDF.mod)\n\n";

}
print "DONE."

