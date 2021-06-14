#!/usr/bin/env perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo
#
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2012 Universidade de São Paulo

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use Storable;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use Bio::SeqIO;
use Bio::Seq;
use Bio::Index::Fasta;
use DBI;

use constant MAX_TRIALS_TO_FETCH=>2;

my ($level,$kuin,$infile,$excludes_line,$outfile,$tmpdir,$db,$fafile);
# Command line named arguments (see function Usage() )
Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "k|k=s"=>\$kuin,
            "i|infile=s"=>\$infile,
            "e|excludes=s"=>\$excludes_line,
            "o|outfile=s"=>\$outfile,
            "t|tmpdir=s"=>\$tmpdir,
            "d|db=s"=>\$db,
            "f|fasta=s"=>\$fafile
    ) or &Usage();


#Log4perl levels
if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}

#DataBase handle (DB connection)
my $dbh;
#Multifasta file index
my $inx;

if ($fafile) {
    
    $LOGGER->logdie("Wrong fasta file ($fafile)") unless (-e $fafile);
    

    unless ( -e $fafile.'.idx' ) {
        # Make an index for one or more fasta files
         
        $inx = Bio::Index::Fasta->new(-filename => $fafile.'.idx',
                                         -write_flag => 1);
        $inx->make_index($fafile);

    } else {
        $inx = Bio::Index::Fasta->new(-filename => $fafile.'.idx');
    }

}

# SQLite database connection if db file is present in command line
# Not used today
if ($db) {
    
    $LOGGER->logdie("Wrong KEGG db file ($db)") unless (-e $db);
    $dbh = DBI->connect("dbi:SQLite:dbname=$db","","");
	
}


# Create a user agent object
use LWP::UserAgent;
my $ua = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");

# Diretory to write files (not so temporary)
$tmpdir||='/tmp';

$LOGGER->logdie("Not found temporary directory ($tmpdir)") unless (-e $tmpdir);

$tmpdir=~s/\/+$//;

# Hash to store lineages for exclusion, i.e. which will be ignored
my %excludes;
if ($excludes_line) {
    chomp($excludes_line);
    foreach my $e (split(/,/, $excludes_line)) {
        $excludes{$e} = undef;
    }
}

# Input IDs
my %kin;


unless ($infile) {
    # Only one ID (using -k)
     
    $LOGGER->logdie("Missing k id or input file with k ids") unless ($kuin);

    $kin{$kuin}->{''} = undef;

    undef($kuin);

} else {
    # Multiple IDs in an input file

    $LOGGER->logdie("Wrong file ($infile)") unless (-e $infile);

    open(IN, "<", $infile) or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        next if ($_=~/^#/);

        # The file can contain ID (obligatory) and a category description
        # The ID can belongs to multiple category, so we used a HASH of HASH references
        my ($id, $desc) = split(/\t/, $_);
        $kin{$id}->{$desc||''} = undef;
    }
    close(IN);

}

# Code to recovery of Kegg Orthology IDs from other IDs, such as pathways (ko), modules (M) or brite (br)
foreach my $id (keys %kin) {
    
    if ($id =~ /^(ko|M)[0-9]{5,5}/) { #ko or M
            
            my ($cat) = $1;

            my $hr_desc = $kin{$id};

            delete($kin{$id});
            
            my $content_k;
            
            # verify if the web content was already saved in a file and that file is not empty
            if ( ( ! -e "$tmpdir/$cat/$id.txt") || ( -z "$tmpdir/$cat/$id.txt") ) {
                # Create a request
                my $req_k = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$id);
                $req_k->content_type('application/x-www-form-urlencoded');
                $req_k->content('query=libwww-perl&mode=dist');
                #sleep(1);

                # Pass request to the user agent and get a response back
                my $res_k = $ua->request($req_k);

                # Check the outcome of the response
                if ($res_k->is_success()) {
                            $content_k = $res_k->content();
                } else {
                    $LOGGER->logwarn( $res_k->status_line );
                }
                
                # save content into a file
                open(CONTENT, ">", "$tmpdir/$cat/$id.txt") or $LOGGER->logdie($!);
                print CONTENT $content_k;
                close(CONTENT);
                
                
            } else {
                # recovery content from a file
                open(CONTENT, "<", "$tmpdir/$cat/$id.txt") or $LOGGER->logdie($!);
                while(<CONTENT>) {
                    $content_k.=$_;
                }
                close(CONTENT);
            }
            
            my $set_k;
            # parse content and extract Kegg Orthology IDs
            foreach my $line_k ( split(/\n/, $content_k) ) {
                
                if (($line_k =~ /^ORTHOLOGY/)||($set_k)) {
                    last if (($set_k)&&($line_k=~/^\S+/));
                    $set_k = 1;
                    my ($ostr) = $line_k=~/^(?:ORTHOLOGY)?\s+(\S+)/;
                    foreach my $oid (split(/,|\+/, $ostr)) {
                        foreach my $desc (keys %{ $hr_desc }) {
                            $kin{$oid}->{$desc} = undef;
                        }
                    }
                }
            }
            
        } elsif ($id =~ /^(br):ko[0-9]{5,5}/) { #br:ko
            
            my $cat = $1;
            
            my $hr_desc = $kin{$id};
            
            delete($kin{$id});
            
            my $content_k;
            
            if ( ( ! -e "$tmpdir/$cat/$id.txt") || ( -z "$tmpdir/$cat/$id.txt") ) {
                # Create a request
                my $req_k = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$id);
                $req_k->content_type('application/x-www-form-urlencoded');
                $req_k->content('query=libwww-perl&mode=dist');
                #sleep(1);
                
                # Pass request to the user agent and get a response back
                my $res_k = $ua->request($req_k);
                
                # Check the outcome of the response
                if ($res_k->is_success()) {
                            $content_k = $res_k->content();
                } else {
                    $LOGGER->logwarn( $res_k->status_line );
                }
                
                open(CONTENT, ">", "$tmpdir/$cat/$id.txt") or $LOGGER->logdie($!);
                print CONTENT $content_k;
                close(CONTENT);
                
            } else {
                open(CONTENT, "<", "$tmpdir/$cat/$id.txt") or $LOGGER->logdie($!);
                while(<CONTENT>) {
                    $content_k.=$_;
                }
                close(CONTENT);       
            }
            
            my $set_k;
            
            foreach my $line_k ( split(/\n/, $content_k) ) {
                
                if ($line_k =~ /^\+D\s+(.*)/) {
                    delete($hr_desc->{''});
                    $hr_desc->{$1} = undef;
                } elsif ($line_k =~ /^[A-Z]\s+(K[0-9]{5,5})\s+\S+.*/) { 
                    my $oid=$1;
                    foreach my $desc (keys %{ $hr_desc }) {
                        $kin{$oid}->{$desc} = undef;
                    }
                }
            }
        
        }
}

my $seqout;

# Object to write the sequences in FASTA format
if ($outfile) {
    $seqout=Bio::SeqIO->new(-file=>'>'.$outfile, 
                            -format=>'FASTA',
                            -flush=>0);
}


my $hr_org;

unless (-e "./keggorganism.dump") {
    
    # Create a request
    my $req = HTTP::Request->new(GET => 'http://rest.kegg.jp/list/organism');
    $req->content_type('application/x-www-form-urlencoded');
    $req->content('query=libwww-perl&mode=dist');
    
    # Pass request to the user agent and get a response back
    my $res = $ua->request($req);
    
    # Check the outcome of the response
    if ($res->is_success()) {
        my $set = undef;
        my $content = $res->content();
        foreach my $line ( split(/\n/, $content) ) {
                my ($id, $code, $name, $lineage) = split(/\t/, $line);
                $hr_org->{$code} = {'id'=>$id,
                                    'name'=>$name,
                                    'lineage'=>$lineage};
        }
        # Save the hash reference into a file
        store $hr_org, "./keggorganism.dump";
    }
    else {
        print $res->status_line, "\n";
    }
    
} else {
    # Retrieve the hash reference from a file
    $hr_org = retrieve("./keggorganism.dump");
}

my %valid_org = %{ $hr_org };
# Exclude (delete) lineages from data structure
# according with -e command line parameter
foreach my $code (keys %{$hr_org}) {
    foreach my $e (keys %excludes) {
        if ($valid_org{$code}) {
            if ($valid_org{$code}->{'lineage'}=~/$e/) {
                delete( $valid_org{$code} );
                last;
            }
        }
    }
}

# Code for testing Hash reference structure and content of
# KEGG Lineage Data 
#foreach my $code (keys %valid_org) {
#    print $code,"\t",$valid_org{$code}->{'name'},"\t",$valid_org{$code}->{'lineage'},"\n";
#}

$LOGGER->logdie("Not found organism list") unless (scalar(keys(%{$hr_org}))>0);

$|=1;
my %OUTPUT;
my %PRINT;

foreach my $kid (keys %kin) {
    
    $kid=~s/\s+//g;

    #   print $kid,"\n";
    my $content_k;
    
    if ($kid !~ /^K[0-9]{5,5}/) {
        
        my $cat = 'uniprot';
        
        print STDERR "Searching $kid in UniProt ... ";
        
        my $seqobj;
        
        if ( ( ! -e "$tmpdir/$cat/$kid.txt") || ( -z "$tmpdir/$cat/$kid.txt") ) {
            
            print STDERR " DOWNLOADING ";
                                
            my $auxseqout = Bio::SeqIO->new( -file   => ">$tmpdir/$cat/$kid.txt", 
                                             -format => 'swiss', 
                                             -flush  => 0 );
            
            # Create a request
            my $req_k = HTTP::Request->new(GET => 'https://www.uniprot.org/uniprot/'.$kid.'.txt');
            $req_k->content_type('application/x-www-form-urlencoded');
            $req_k->content('query=libwww-perl&mode=dist');
            #sleep(1);
            
            # Pass request to the user agent and get a response back
            my $res_k = $ua->request($req_k);
            
            # Check the outcome of the response
            if ($res_k->is_success()) {
                $content_k = $res_k->content();
                
                my $auxseqin = Bio::SeqIO->new( -string=> $content_k, 
                                                -format => 'swissdriver');
                
                $seqobj = $auxseqin->next_seq();
                
                $auxseqout->write_seq($seqobj);

            } else {
                $LOGGER->logwarn( $res_k->status_line );
                next;
            }
            

        } else {
            
                my $auxseqin = Bio::SeqIO->new( -file   => "$tmpdir/$cat/$kid.txt", 
                                                -format => 'swissdriver' );
                $seqobj = $auxseqin->next_seq();
        }
        
        my @n;
        my $name=$kid;
        
        #my @gn;
        #while ($content_k=~/^GN\s+(.*)/mg) {
        #    push(@gn, $1);
        #}
        #my $joingn=join(' ',@gn);
        #if ($joingn) {
        #    my @n;
        #    while($joingn=~/\S+=([^;]+)/g) {
        #        my $eachn = $1;
        #        $eachn=~s/\s\{[^\}]+\}//;
        #        push(@n, $eachn);
        #    }
        #    $name=join(', ', @n);
        #}
        #
        #my $seq;
        #do {
        #    local *STDERR;
        #    open(STDERR, ">", "./getInfoByK.log.err.txt") or $LOGGER->logdie($!);
        #    my $seqin = Bio::SeqIO->new(-string => $content_k, -format => 'swissdriver', -verbose=>-1);
        #    $seq = $seqin->next_seq();
        #};
        
        my @species = $seqobj->species()->classification();
        
        my $lineage = join(";", reverse(@species));
        
        my $set = 1;
        foreach my $e (keys %excludes) {
                if ($lineage=~/$e/) {
                    undef($set);
                    last;
                }
        }
        
        if ($set) {
            print STDERR " FOUND\n";
            
            my $definition = $seqobj->description();
            $definition=~s/;$//;
            
            if ($kin{$kid}) {
                foreach my $d (keys %{ $kin{$kid} }) {
                    $PRINT{$kid} = ''.(($d) ? $d."\t" : '').$kid."\t".$name."\t".$definition;
                }                    
            }
            
            if ($seqout) {
                $seqobj->display_id($kid);
                $OUTPUT{$kid}->{$kid} = $seqobj;
            }
        }
        
    } else {
        
        my $cat = 'KO';

        my $trial_ok= undef;
        my $trial=1;

        while ( (! defined $trial_ok) && ($trial <= MAX_TRIALS_TO_FETCH) ) {
            if ($trial>1) {
                print STDERR "[Trial number: $trial]\n";
            }
            $trial_ok=1;
            
            if ( ( ! -e "$tmpdir/$cat/$kid.txt") || ( -z "$tmpdir/$cat/$kid.txt") ) {
                
                # Create a request
                my $req_k = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$kid);
                $req_k->content_type('application/x-www-form-urlencoded');
                $req_k->content('query=libwww-perl&mode=dist');
                #sleep(1);
                
                # Pass request to the user agent and get a response back
                my $res_k = $ua->request($req_k);
                
                # Check the outcome of the response
                if ($res_k->is_success()) {
                        $content_k = $res_k->content();

                        open(CONTENT, ">", "$tmpdir/$cat/$kid.txt") or $LOGGER->logdie($!);
                        print CONTENT $content_k;
                        close(CONTENT);

                } else {
                    $LOGGER->logdie( $kid.":".$res_k->status_line );
                    $trial++ if ($trial_ok);
                    $trial_ok=undef;
                    next;
                }
                
            } else {
                open(CONTENT, "<", "$tmpdir/$cat/$kid.txt") or $LOGGER->logdie($!);
                while(<CONTENT>) {
                    $content_k.=$_;
                }
                close(CONTENT);
            }
            
            my $set_k = undef;
            my $name_k = undef;
            my $definition_k = undef;
            my $name;
            my $definition;
            my $set_print=undef;
            
            foreach my $line_k ( split(/\n/, $content_k) ) {
                if (($line_k =~ /^NAME/)||($name_k)) {
                    if (($name_k)&&($line_k=~/^\S+/)) {
                        $name_k=undef;
                    } else {
                        $name_k=1;
                        my ($tmpname) = $line_k=~/^(?:NAME)?\s+(\S+.*)/;
                        if ($name) {
                            $name.=' ';
                        }
                        $name.=$tmpname;
                    }                
                }

                if (($line_k =~ /^DEFINITION/)||($definition_k)) {
                    if (($definition_k)&&($line_k=~/^\S+/)) {
                        $definition_k=undef;
                    } else {
                        $definition_k=1;
                        my ($tmpdef) = $line_k=~/^(?:DEFINITION)?\s+(\S+.*)/;
                        if ($definition) {
                            $definition.=' ';
                        }
                        $definition.=$tmpdef;
                    } 
                }
                
                if (($line_k =~ /^GENES/)||($set_k)) {
                    last if (($set_k)&&($line_k=~/^\S+/));
                    $set_k = 1;
                    my ($db, $acc_list) = $line_k=~/^(?:GENES)?\s+(\S+):\s+(.+)/;
                    if ( (exists $valid_org{lc($db)}) || (! exists $hr_org->{lc($db)}) ) {

                        if (! exists $hr_org->{lc($db)}) {
                            $LOGGER->logwarn("Not found organism code \"".lc($db)."\"");
                        }

                        unless (defined $set_print) {

                            if ($kin{$kid}) {
                                foreach my $d (keys %{ $kin{$kid} }) {
                                    $PRINT{$kid} = ''.(($d) ? $d."\t" : '').$kid."\t".$name."\t".$definition;
                                }                                            
                            }

                            $set_print=1;
                        }
                        
                        my $seqobj;
                        my $seqid;
                        
                        if ($seqout) {

                            my $new_acc_list=""; 
                            my @parenthesis;
                            my @aux=split(//, $acc_list);

                            for my $j (@aux) {  
                                if($j eq "(" ) { 
                                    push(@parenthesis,1); 
                                } elsif ($j eq ")") { 
                                    pop(@parenthesis); 
                                } else { 
                                    $new_acc_list.=$j if (scalar(@parenthesis)==0);  
                                }   
                            }

                            foreach my $acc ( split(/ /, $new_acc_list) ) {
                                
                                #$acc=~s/\([^\)]+\)+//g;
                                $acc=~s/^\s+//;
                                $acc=~s/\s+$//;
                                $LOGGER->logdie("Wrong Accession: $acc") if ($acc!~/^[A-Za-z0-9\_\-\.]+$/);
                                $seqid=lc($db).':'.$acc;
                                
                                die "$kid/$seqid" if ($acc=~/[\(\)]/);   
                                next if (exists $OUTPUT{$kid}->{$seqid});

                                my $cat = lc($db);

                                if ( ( ! -e "$tmpdir/$cat/$acc.fa") || ( -z "$tmpdir/$cat/$acc.fa") ) {
                                    
                                    mkdir("$tmpdir/$cat/") unless (-d "$tmpdir/$cat");
                                
                                    my $auxseqout = Bio::SeqIO->new( -file   => ">$tmpdir/$cat/$acc.fa", 
                                                                     -format => 'FASTA', 
                                                                     -flush  => 0 );
                                    
                                    if (defined $inx) {
                                        $seqobj = $inx->fetch($seqid);
                                    }

                                    unless ($seqobj) {
                                        my $req_s = HTTP::Request->new(GET => 'http://rest.kegg.jp/get/'.$seqid.'/aaseq');
                                        $req_s->content_type('application/x-www-form-urlencoded');
                                        $req_s->content('query=libwww-perl&mode=dist');
                                        #sleep(1);
                                        
                                        my $res_s= $ua->request($req_s);
                                        my $content_s;
                                        
                                        if ($res_s->is_success()) {
                                            $content_s = $res_s->content();

                                            my $auxseqin = Bio::SeqIO->new(-string=> $content_s, -format => 'FASTA');
                                        
                                            $seqobj = $auxseqin->next_seq();

                                        } else{
                                            $LOGGER->logwarn( "Trying to get [$kid] ".$seqid."  (". $res_s->status_line.")" );
                                            $trial++ if ($trial_ok);
                                            $trial_ok=undef;
                                            next;
                                        }
                                        
                                        
                                    } else {
                                        $seqobj->desc("$kid $definition | Recovered by searching KOBAS sequences");
                                    }
                                    
                                    $auxseqout->write_seq($seqobj);
                                
                                } else {
                                    
                                    my $auxseqin = Bio::SeqIO->new( -file   => "$tmpdir/$cat/$acc.fa", 
                                                                    -format => 'FASTA' );
                                    
                                    $seqobj = $auxseqin->next_seq();
                                    
                                }
                                
                                $seqobj->display_id($kid.':'.$seqid);
                                $seqobj->desc("");
                                 
                                $OUTPUT{$kid}->{$seqid} = $seqobj if ($seqobj);
                            }
                        }
                    }
                }
            }
            if (! defined $trial_ok) {
                unlink("$tmpdir/$cat/$kid.txt");
            }                
        }

    }        
}
    
foreach my $kid (keys %PRINT) {
        print $PRINT{$kid},"\n";
}

if ($seqout) {
    foreach my $kid (keys %OUTPUT) {
        foreach my $seqid (keys %{$OUTPUT{$kid}}) {
            $seqout->write_seq( $OUTPUT{$kid}->{$seqid} );
        }
    }
    $seqout = undef;
}

if ($dbh) {
    $dbh->disconnect();
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -k      --k         KEGG Orthology (KO), KEGG reference pathway (ko), KEGG module (M) or UniProt ID
        -i      --in        Input file where: the first column has KEGG Orthology (KO), KEGG reference pathway (ko), KEGG module (M) or UniProtID 
                                              the second column (optional) has description
        -e      --excludes  Exclude lineages (comma separated names)
        -o      --outfile   Output FASTA file
        -t      --tmpdir    Temporary directory [Default: /tmp]
		-d		--db		KOBAS db SQLite file (ko.db)
        -f      --fasta     KOBAS FASTA file (ko.pep.fasta)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

