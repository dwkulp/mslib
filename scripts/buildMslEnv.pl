#!/usr/bin/perl
use File::Fetch;
use File::Basename;

use Getopt::Long;
&GetOptions(
	"name=s"        => \$name,
        "debug"         => \$debug,
        "j=s"           => \$J
    );

$J=1 if (!$J);
if (!$name){
    if ($#ARGV != -1){
	$name = $ARGV[0];
    } else {
	$name = "extlib";
    }

}
$pwd=`pwd`; $pwd =~ s/\n//;
print "External library directory is: ".$pwd."/".$name."_libs/\n";

my %libs;
open(TMP, "libs.txt") or die "Couldn't open libs.txt: $?\n";
while ($line = <TMP>){
    $line =~ s/\n//;
    next if (substr(0,1) eq "#");

    if ($line =~ /^(\S+):/){
	
	$lib = $1;
	$version = "";
	if ($line =~ /^(\S+):\s*(\S+)$/){
	    $version = $2;
	}

	$libs{$lib} = $version;
    }
}
close(TMP);

# bash
#export MSL_GSL=T
#export MSL_GLPK=T
#export MSL_BOOST=T
#export MSL_OPENMP=T
#export MSL_EXTERNAL_LIB_DIR=/usr/lib
#export MSL_EXTERNAL_INCLUDE_DIR=/usr/include

# tcsh 
# setenv MSL_GSL T
# setenv MSL_GLPK T
# setenv MSL_BOOST T
# setenv MSL_R T
# setenv MSL_OPENMP T
# setenv MSL_EXTERNAL_LIB_DIR /usr/lib
# setenv MSL_EXTERNAL_INCLUDE_DIR /usr/include

# Make a lib_src directory
system("mkdir -p $name"."_src");
system("mkdir -p $name"."_libs");
system("mkdir -p $name"."_incs $name"."_incs/gsl");
system("cp libs.txt $name"."_libs/");
system("cp libs.txt $name"."_incs/");

my $pwd=`pwd`; $pwd =~ s/\n//;
my $bash_lines ="export MSL_EXTERNAL_LIB_DIR=$pwd/$name"."_libs\n";
$bash_lines .= "export MSL_EXTERNAL_INCLUDE_DIR=$pwd/$name"."_incs\n";
my $tcsh_lines ="setenv MSL_EXTERNAL_LIB_DIR $pwd/$name"."_libs\n";
$tcsh_lines .= "setenv MSL_EXTERNAL_INCLUDE_DIR $pwd/$name"."_incs\n";

# for each keys
foreach $lib (keys %libs){
    if ($libs{$lib} eq ""){
	print "Not building library $lib\n";
	$bash_lines .= "export ".$lib."=F\n";
	$tcsh_lines .= "setenv ".$lib." F\n";
	next;
    } 

    $bash_lines .= "export ".$lib."=T\n";
    $tcsh_lines .= "setenv ".$lib." T\n";

    if (!$debug && $lib eq "MSL_GSL"){
	my $ff = File::Fetch->new(uri => $libs{$lib});
	my $where = $ff->fetch(to => "./$name"."_src/") or die $ff->error;
	my $fname = $ff->file;
	my $base = "";
	if ($fname =~ /(\S+).tar.gz/){
	    $base = $1;
	} else {
	    print "ERROR MSL_GSL building tar.gz not name $libs{$lib},$fname,$base\n";
	    exit(1);
	}
	system("cd ./$name"."_src; tar xzvf $fname");
	system("cd ./$name"."_src/$base; ./configure --prefix=`pwd`/../../$name"."_libs/");
	system("cd ./$name"."_src/$base; make -j$J");
	system("cd ./$name"."_src/$base; cp .libs/libgsl.a ../../$name"."_libs/");
	system("cd ./$name"."_src/$base; cp cblas/.libs/libgslcblas.a ../../$name"."_libs/");
	system("cd ./$name"."_src/$base; cp -H gsl/* ../../$name"."_incs/gsl");
    }

    if (!$debug && $lib eq "MSL_BOOST"){

	# Instal bzp2 lib for boost::iostreams
	#system("sudo apt-get install libbz2-dev");

	$fname = fileparse($libs{$lib});
	my $ff = File::Fetch->new(uri => $libs{$lib});
	print "Fetching $libs{$lib}, $name\n";
	print "Directory exists: ".(-d "./$name"."_src/")."\n";
	my $where = $ff->fetch(to => "./$name"."_src/") or die $ff->error;
	my $fname = $ff->file;
        my $base = "";
	if ($fname =~ /(\S+).tar.gz/){
	    $base = $1;
	} else {
  	    print "ERROR MSL_BOOST building tar.gz not name $libs{$lib},$fname,$base\n";
	    exit(1);
	}
	system("cd ./$name"."_src; tar xzvf $fname");


#	system("cd ./$name"."_src/; svn co $libs{$lib}");
	system("cd ./$name"."_src/$base; ./bootstrap.sh --prefix=`pwd`/../../$name"."_libs/ --with-libraries=iostreams,serialization,regex,filesystem");
	system("cd ./$name"."_src/$base; ./b2");
	system("cd ./$name"."_src/$base; cp ./stage/lib/libboost_*.a ../../$name"."_libs/");
	system("cd ./$name"."_src/$base; cp -r boost ../../$name"."_incs/");
    }

    if (!$debug && $lib eq "MSL_GLPK"){
	my $ff = File::Fetch->new(uri => $libs{$lib});
	my $where = $ff->fetch(to => "./$name"."_src/") or die $ff->error;
	my $fname = $ff->file;
	my $base = "";
	if ($fname =~ /(\S+).tar.gz/){
	    $base = $1;
	} else {
	    print "ERROR MSL_GLPK building tar.gz not name $libs{$lib},$fname,$base\n";
	    exit(1);
	}

	system("cd ./$name"."_src; tar xzvf $fname");
	system("cd ./$name"."_src/$base; ./configure --prefix=`pwd`/dir_install/ --disable-dl --disable-odbc --disable-mysql --enable-static --without-gmp --without-zlib; echo `pwd`/dir_install");
	
	system("cd ./$name"."_src/$base; make -j$J");
	system("cd ./$name"."_src/$base; make install");
	system("cd ./$name"."_src/$base; cp dir_install/lib/libglpk.a ../../$name"."_libs/");
	system("cd ./$name"."_src/$base; cp dir_install/include/glpk.h ../../$name"."_incs/");


    }

    

    
}

print <<TXT;
############# FOR BASH .bashrc/.bash_profile  ################
$bash_lines

############# FOR TCSH .tcshrc  ################
$tcsh_lines

TXT
