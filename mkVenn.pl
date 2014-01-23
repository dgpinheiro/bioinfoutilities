#!/usr/bin/perl
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

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, %infiledef, %incolordef, $tifffile, $outfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\%infiledef,
            "c|incolor=s"=>\%incolordef,
            "t|tifffile=s"=>\$tifffile,
            "o|outfile=s"=>\$outfile
    ) or &Usage();

$LOGGER->logdie("Missing input file") unless (scalar keys %infiledef);

$LOGGER->logdie("Missing output file") unless ($outfile);

my $outfile_up = $outfile;
$outfile_up=~s/(\.[^\.]+)$/_up$1/;

my $outfile_down = $outfile;
$outfile_down=~s/(\.[^\.]+)$/_down$1/;

my $outfile_ctreg = $outfile;
$outfile_ctreg=~s/(\.[^\.]+)$/_ctreg$1/;

foreach my $def (keys %infiledef) {
    $LOGGER->logdie("Wrong input file ($infiledef{$def})") unless (-e $infiledef{$def});
}

$LOGGER->logdie("Missing color(s)") unless (scalar keys %incolordef);

foreach my $k (keys %infiledef) {
    $LOGGER->logdie("Please don't use space for definitions ($k)") if ($k=~/\s/);
}

$LOGGER->logdie("Missing .tiff file") unless ($tifffile);

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

use File::Basename;

use R;
use Cwd;

&R::initR("--silent");
&R::library("RSPerl");

&R::eval('library(package="VennDiagram")');
&R::eval('library(package="colorspace")');
&R::eval('library(package="gridExtra")');
&R::eval('lst.up <<- list()');
&R::eval('lst.down <<- list()');
&R::eval('lst.all <<- list()');
my @alldefs = sort keys %infiledef;
my @allcolors = @incolordef{ @alldefs };

&R::eval("fillcolor <<- data.frame(ID=c(".join(',', map { '"'.$_.'"' } @alldefs ).'), Color=c('.join(',', map { '"'.$_.'"' } @allcolors ).'))
          rownames(fillcolor) <<- fillcolor$ID
');
&R::eval('catcolor <<- fillcolor
          rownames(catcolor) <<- catcolor$ID
          catcolor[,"Color"] <<- apply(catcolor, 1, function(x) { 
                                                        hex(mixcolor(0.35, hex2RGB(
                                                        rgb(red=col2rgb(x["Color"])["red",],
                                                            green=col2rgb(x["Color"])["green",],
                                                            blue=col2rgb(x["Color"])["blue",],
                                                            maxColorValue=255
                                                        )),
                                                        RGB(0, 0, 0))) 
      
            })');

#&R::eval('save.image(file = "/tmp/TESTEBegin.RData", version = NULL, safe = TRUE)');



my %analysis;
foreach my $def (@alldefs) {

    $LOGGER->logdie("Not found color definition for $def") unless ($incolordef{$def});
    
    my $infile = $infiledef{$def};

    my $basename = basename($infile);
    $basename=~s/\..*$//;
    
    print "Processing ($def): $basename ...\n";

    my @data_up;
    my @data_down;
    my @data;

    open(IN, "<", $infile) or $LOGGER->logdie($!);

    while(<IN>) {
        chomp;
        my ($item, $fold) = split(/\t/, $_);
        if ($fold > 0) {
            push(@data_up, $item);
        }
        else {
            push(@data_down, $item);
        }
        push(@data, $item);
    }

    close(IN);
    
    &R::eval("lst.up[['$def']] <<- c(".join(',', map { '"'.$_.'"' } @data_up).')');

    &R::eval("lst.down[['$def']] <<- c(".join(',', map { '"'.$_.'"' } @data_down).')');

    &R::eval("lst.all[['$def']] <<- c(".join(',', map { '"'.$_.'"' } @data).')');
   

    $analysis{$basename} = undef;
}

&R::eval('
    get.test.cat <<- function(r, tlst) {
        itm <- r["Item"]
        unlst <- unlist(lapply(tlst, function(r,n) { return(any(as.character(r)==n)) }, itm))
        tmp<-c()
        tmp[unlst] <- 1
        tmp[!unlst] <- 0
        return(paste(tmp,collapse=";"))
    }
    
    p1 <<- NULL
    p2 <<- NULL
    p3 <<- NULL
        
    item.up <<- unique(unlist(lst.up))
    item.down <<- unique(unlist(lst.down))
    
    item.all <<- unique(unlist(lst.all))

    if (length( item.all ) > 0) {

        test.df.all <<- data.frame(Item=item.all,stringsAsFactors=FALSE)

        test.df.all[[paste(names(lst.all),collapse=";")]] <<- apply(test.df.all, 1, get.test.cat, lst.all)
            
        if (length(item.up) > 0) {
    
            test.df.up <<- test.df.all[ apply(test.df.all, 1,
                                        function(x,t.up,t.down) { return( ((any(t.up==as.character(x["Item"])))&(!any(t.down==as.character(x["Item"])))) ) } 
                                        , item.up, item.down), ]


            write.table(test.df.up,file="'.$outfile_up.'",row.names=FALSE,quote=FALSE,sep="\t")

            header = as.character(unlist(strsplit(colnames(test.df.up)[2],";")))
            lst.up <- list()
            rownames(test.df.up) <- test.df.up$Item
            for (i in rownames(test.df.up) ) {
                tstbool <- unlist(strsplit(test.df.up[i,2],";")) == "1"
                names(tstbool) <- header
                for (h in header) {
                    if (tstbool[h]) {
                        lst.up[[h]][length( lst.up[[h]] )+1] <- i
                    }
                }
            }

            up.grob.list <<- venn.diagram(x = lst.up, filename = NULL, 
                                         col = "transparent", fill = as.character(fillcolor[names(lst.up), "Color"]), cat.col = as.character(catcolor[names(lst.up), "Color"]), alpha = 0.50, 
                                         cex = 3, fontfamily = "serif", fontface = "bold",
                                         cat.cex = 1.5, cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", margin = 0.1, main="Up-regulated") 

            p1 <<- grid.grabExpr(grid.draw(up.grob.list))
        } else {
            p1 <<- grid.grabExpr(grid.draw( linesGrob() ))
        }

        if (length(item.down) > 0) {
            
            test.df.down <<- test.df.all[ apply(test.df.all, 1,
                                        function(x,t.up,t.down) { return( ((!any(t.up==as.character(x["Item"])))&(any(t.down==as.character(x["Item"])))) ) } 
                                        , item.up, item.down), ]

        
            write.table(test.df.down,file="'.$outfile_down.'",row.names=FALSE,quote=FALSE,sep="\t")

            header = as.character(unlist(strsplit(colnames(test.df.down)[2],";")))
            lst.down <- list()
            rownames(test.df.down) <- test.df.down$Item
            for (i in rownames(test.df.down) ) {
                tstbool <- unlist(strsplit(test.df.down[i,2],";")) == "1"
                names(tstbool) <- header
                for (h in header) {
                    if (tstbool[h]) {
                        lst.down[[h]][length( lst.down[[h]] )+1] <- i
                    }
                }
            }

            down.grob.list <<- venn.diagram(x = lst.down, filename = NULL, 
                                         col = "transparent", fill = as.character(fillcolor[names(lst.down), "Color"]), cat.col = as.character(catcolor[names(lst.down), "Color"]), alpha = 0.50, 
                                         cex = 3, fontfamily = "serif", fontface = "bold",
                                         cat.cex = 1.5, cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", margin = 0.1, main="Down-regulated") 

            p2 <<- grid.grabExpr(grid.draw(down.grob.list))
                                                        
        } else {
            p2 <<- grid.grabExpr(grid.draw( linesGrob() ))
        }

        if ( (length(item.down) > 0)&(length(item.up) > 0) ) {

            test.df.ctreg <<- test.df.all[ apply(test.df.all, 1,
                                        function(x,t.up,t.down) { return( ((any(t.up==as.character(x["Item"])))&(any(t.down==as.character(x["Item"])))) ) } 
                                        , item.up, item.down), ]
            if (dim(test.df.ctreg)[1]>0) {
                write.table(test.df.ctreg,file="'.$outfile_ctreg.'",row.names=FALSE,quote=FALSE,sep="\t")
           
                header = as.character(unlist(strsplit(colnames(test.df.ctreg)[2],";")))
                lst.ctreg <- list()
                rownames(test.df.ctreg) <- test.df.ctreg$Item
                for (i in rownames(test.df.ctreg) ) {
                    tstbool <- unlist(strsplit(test.df.ctreg[i,2],";")) == "1"
                    names(tstbool) <- header
                    for (h in header) {
                        if (tstbool[h]) {
                            lst.ctreg[[h]][length( lst.ctreg[[h]] )+1] <- i
                        }
                    }
                }
         
                ctreg.grob.list <<- venn.diagram(x = lst.ctreg, filename = NULL, 
                                             col = "transparent", fill = as.character(fillcolor[names(lst.ctreg), "Color"]), cat.col = as.character(catcolor[names(lst.ctreg), "Color"]), alpha = 0.50, 
                                             cex = 3, fontfamily = "serif", fontface = "bold",
                                             cat.cex = 1.5, cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", margin = 0.1, main="Contra-regulated", ext.text=TRUE, cat.default.pos="outer", offset=0.5) 
                
                p3 <<- grid.grabExpr(grid.draw(ctreg.grob.list))
            } else {
                p3 <<- grid.grabExpr(grid.draw( linesGrob() ))
            }
        } else {
            p3 <<- grid.grabExpr(grid.draw( linesGrob() ))
        }
        
    }

    #save.image(file = "/tmp/mkvenn_debug.RData", version = NULL, safe = TRUE)

    if ( (!is.null(p1)) && (!is.null(p2)) ) {
        tiff(filename = "'.$tifffile.'", height = 3000, width = 9000, compression="lzw", res=300, units="px", type="cairo")
        grid.arrange( p1, p3, p2, 
            ncol = 3,
            main = textGrob( "Venn Diagram of Differentially Expressed Genes", vjust = 1, gp = gpar( fontface = "bold", cex = 1.5 )),
            widths  = unit( 1 / 3, "npc" ),
            heights = unit( 0.8,   "npc" )
        )
        dev.off()
    }




');

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help       Help
        -l      --level      Log level [Default: FATAL]
        -i      --infile     Input definition with .txt file processed by genLists.pl
        -c      --incolor    Input definition with color
        -t      --tifffile   .tiff file
        -o      --outfile    Output file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

