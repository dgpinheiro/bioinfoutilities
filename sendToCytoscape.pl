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

use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $name, $outfile, $pdffile, $validfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "n|name=s"=>\$name,
            "o|outfile=s"=>\$outfile,
            "j|pdffile=s"=>\$pdffile,
            "v|validfile=s"=>\$validfile
    ) or &Usage();


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

$LOGGER->logdie("Missing input file") unless ($infile);
$LOGGER->logdie("Wrong input file ($infile)") unless (-e $infile);

$name||=basename($infile);


use R;
use Cwd;

&R::initR("--silent");
&R::library("RSPerl");


&R::eval('library(package="igraph")');
&R::eval('library(package="plyr")');
&R::eval('library(package="RCytoscape")');

&R::eval('
sendToCytoscape <<- function(networkName="network name", graphFile, validFile) {
    
    # Read a data set. 
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet <- read.table(graphFile, header = FALSE, sep = "\t")
    validSet <- read.table(validFile, header = FALSE, sep = "\t")

    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD <- simplify(graph.data.frame(dataSet, directed=FALSE))
    
    # Print number of nodes and edges
    # vcount(gD)
    # ecount(gD)
    
    ############################################################################################
    # Calculate some node properties and node similarities that will be used to illustrate 
    # different plotting abilities
    
    # Calculate igraph::degree for all nodes
    degAll <- igraph::degree(gD, v = V(gD), mode = "all")
    
    # Calculate betweenness for all nodes
    betAll <- betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
    betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
    rm(betAll)
    
    #Calculate Dice similarities between all pairs of nodes
    dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
    
    ############################################################################################
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    findType <- function(name) {
        if (grepl("^[a-z]{3}-(let|bantam|miR)", name)) {
            return("miRNA");
        } else {
            return("mRNA");
        }
    }
    typeAll <- unlist(lapply(V(gD)$name, findType))
    names(typeAll) <-  V(gD)$name
    
    findStatus <- function(name, validset) {
        if (name %in% validset) {
            return("Y")
        } else {
            return("N")
        }
    }
    statusAll <- unlist(lapply(V(gD)$name, findStatus, levels(validSet$V1)))
    names(statusAll) <-  V(gD)$name

    gD <- set.vertex.attribute(gD, "degree", index = V(gD), value = degAll)
    gD <- set.vertex.attribute(gD, "betweenness", index = V(gD), value = betAll.norm)
    gD <- set.vertex.attribute(gD, "moleculeType", index = V(gD), value = typeAll)
    gD <- set.vertex.attribute(gD, "Status", index = V(gD), value = statusAll)
    # Check the attributes
    # summary(gD)
    
    F1 <- function(x) {data.frame(V4 = dsAll[which(V(gD)$name == as.character(x$V1)), which(V(gD)$name == as.character(x$V2))])}
    dataSet.ext <- ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))
    
    gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0)
    gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
    E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
    
    # Check the attributes
    # summary(gD)
    
    ####################################
    # Print network in Cytoscape
    # This requires RCytoscape package and CytoscapeRPC plugin
    
    
    gD.cyt <- igraph.to.graphNEL(gD)
    
    # We have to create attributes for graphNEL
    # We\'ll keep the same name, so the values are passed from igraph
    
    gD.cyt <- initNodeAttribute(gD.cyt, "degree", "numeric", 0)
    gD.cyt <- initNodeAttribute(gD.cyt, "moleculeType", "char", "") 
    gD.cyt <- initNodeAttribute(gD.cyt, "Status", "char", "N") 
    gD.cyt <- initNodeAttribute(gD.cyt, "betweenness", "numeric", 0) 
    gD.cyt <- initEdgeAttribute(gD.cyt, "weight", "integer", 0)
    gD.cyt <- initEdgeAttribute(gD.cyt, "similarity", "numeric", 0)
    
    # Now we can create a new graph window in cytoscape
    # Be sure that CytoscapeRPC plugin is activated
    gDCW <- new.CytoscapeWindow(networkName, graph = gD.cyt, overwriteWindow = TRUE)
    
    # We can display graph, with defaults color/size scheme
    displayGraph(gDCW)
    
    # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
    cy <- CytoscapeConnection()
    hlp <-getLayoutNames(cy)
    # We\'ll select the layour number 18 - "fruchterman-rheingold" layout 
    # See properties for the given layout
    # getLayoutPropertyNames(cy, hlp[18])
    # Apply values to some of the properties and plot the layout
    setLayoutProperties (gDCW, hlp[18], list (edge_attribute = "similarity", iterations = 1000))
    
    layoutNetwork(gDCW, hlp[18])
    
    # Figure 1 made here
    
    # Now, we can define our own default color/size scheme
    setDefaultBackgroundColor(gDCW, "#FFFFFF")
    setDefaultEdgeColor(gDCW, "#CDC9C9")
    setDefaultEdgeLineWidth(gDCW, 4)
    setDefaultNodeBorderColor(gDCW, "#000000")
    setDefaultNodeBorderWidth(gDCW, 3)
    #setDefaultNodeShape(gDCW, "ellipse")
    setNodeShapeRule (gDCW, node.attribute.name="moleculeType",
            c("miRNA", "mRNA"), c("ellipse", "rect"))
    setNodeBorderColorRule(gDCW, node.attribute.name="Status", c("Y", "N"),   c("#00FF00", "#000000"), mode="lookup", default.color="#000000")
    setNodeBorderWidthRule (gDCW, node.attribute.name="Status", c("Y", "N"), c(24, 3))
    setDefaultNodeColor(gDCW, "#87CEFA")

    setDefaultNodeSize( gDCW, 80 ) 
    setDefaultNodeFontSize( gDCW, 50 )

    setDefaultNodeLabelColor( gDCW, "#000000" )
    
    # And we can replot it 
    redraw(gDCW)       
    
    # Figure 3 made here
    
    # Finally, we can define rules for node colors, node sizes, and edge colors
    setNodeColorRule(gDCW, "degree", c(min(degAll), mean(degAll), max(degAll)), c("#F5DEB3", "#FFA500", "#FF7F50", "#FF4500", "#FF0000"), mode = "interpolate")

    setNodeSizeRule(gDCW, "betweenness", c(min(betAll.norm), mean(betAll.norm), max(betAll.norm)), c(50, 65, 75, 85, 100), mode = "interpolate")

    setEdgeColorRule(gDCW, "weight", c(min(as.numeric(dataSet.ext$V3)), mean(as.numeric(dataSet.ext$V3)), max(as.numeric(dataSet.ext$V3))), c("#CCFFFF", "#00FFFF", "#00FF7F", "#228B22", "#006400"), mode="interpolate")
    
    for (node in names(statusAll[which(statusAll=="N")])) { 
        setNodeLabelOpacityDirect (gDCW, node, 90) 
    } 

    redraw (gDCW)

    return(gDCW)
}
');

&R::eval('
# Clear workspace 
rm(list = ls())
# Send
gDCW <<- sendToCytoscape(networkName="'.$name.'", graphFile="'.$infile.'", validFile="'.$validfile.'")
');


if ($outfile) {
    &R::eval('tbl.noa <<- getAllNodeAttributes (gDCW, onlySelectedNodes=FALSE)
    write.table(
        rbind(
            subset(tbl.noa[order(tbl.noa$degree, decreasing=TRUE),], moleculeType=="miRNA"),
            subset(tbl.noa[order(tbl.noa$degree, decreasing=TRUE),], moleculeType=="mRNA")   
        )
        ,
        "'.$outfile.'",row.names=FALSE,quote=FALSE,sep="\t")
    ');
}

if ($pdffile) {
     &R::eval('saveImage (gDCW, "'.$pdffile.'", "pdf", 1)');
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
        -i      --infile    Input file
        -n      --name      Network name
        -o      --outfile   Text file with node' attributes, such as degree
        -p      --pdffile   PDF file
        -v      --validfile Validated nodes file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

