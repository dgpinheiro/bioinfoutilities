<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html lang="en" dir="ltr">
<head>
<title>TrimBWAstyle.pl - Bioinformatics Core Wiki</title>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<meta name="generator" content="MediaWiki 1.16.1" />
<link rel="shortcut icon" href="/favicon.ico" />
<link rel="search" type="application/opensearchdescription+xml" href="/opensearch_desc.php" title="Bioinformatics Core Wiki (en)" />
<link rel="alternate" type="application/atom+xml" title="Bioinformatics Core Wiki Atom feed" href="/index.php?title=Special:RecentChanges&amp;feed=atom" />
<link rel="stylesheet" href="/skins/common/shared.css?270" media="screen" />
<link rel="stylesheet" href="/skins/common/commonPrint.css?270" media="print" />
<link rel="stylesheet" href="/skins/monobook/main.css?270" media="screen" />
<!--[if lt IE 5.5000]><link rel="stylesheet" href="/skins/monobook/IE50Fixes.css?270" media="screen" /><![endif]-->
<!--[if IE 5.5000]><link rel="stylesheet" href="/skins/monobook/IE55Fixes.css?270" media="screen" /><![endif]-->
<!--[if IE 6]><link rel="stylesheet" href="/skins/monobook/IE60Fixes.css?270" media="screen" /><![endif]-->
<!--[if IE 7]><link rel="stylesheet" href="/skins/monobook/IE70Fixes.css?270" media="screen" /><![endif]-->
<link rel="stylesheet" href="/index.php?title=MediaWiki:Common.css&amp;usemsgcache=yes&amp;ctype=text%2Fcss&amp;smaxage=18000&amp;action=raw&amp;maxage=18000" />
<link rel="stylesheet" href="/index.php?title=MediaWiki:Print.css&amp;usemsgcache=yes&amp;ctype=text%2Fcss&amp;smaxage=18000&amp;action=raw&amp;maxage=18000" media="print" />
<link rel="stylesheet" href="/index.php?title=MediaWiki:Monobook.css&amp;usemsgcache=yes&amp;ctype=text%2Fcss&amp;smaxage=18000&amp;action=raw&amp;maxage=18000" />
<link rel="stylesheet" href="/index.php?title=-&amp;action=raw&amp;maxage=18000&amp;gen=css" />
<script>
var skin="monobook",
stylepath="/skins",
wgUrlProtocols="http\\:\\/\\/|https\\:\\/\\/|ftp\\:\\/\\/|irc\\:\\/\\/|gopher\\:\\/\\/|telnet\\:\\/\\/|nntp\\:\\/\\/|worldwind\\:\\/\\/|mailto\\:|news\\:|svn\\:\\/\\/",
wgArticlePath="/index.php/$1",
wgScriptPath="",
wgScriptExtension=".php",
wgScript="/index.php",
wgVariantArticlePath=false,
wgActionPaths={},
wgServer="http://wiki.bioinformatics.ucdavis.edu",
wgCanonicalNamespace="",
wgCanonicalSpecialPageName=false,
wgNamespaceNumber=0,
wgPageName="TrimBWAstyle.pl",
wgTitle="TrimBWAstyle.pl",
wgAction="view",
wgArticleId=321,
wgIsArticle=true,
wgUserName=null,
wgUserGroups=null,
wgUserLanguage="en",
wgContentLanguage="en",
wgBreakFrames=false,
wgCurRevisionId=2150,
wgVersion="1.16.1",
wgEnableAPI=true,
wgEnableWriteAPI=true,
wgSeparatorTransformTable=["", ""],
wgDigitTransformTable=["", ""],
wgMainPageTitle="Main Page",
wgFormattedNamespaces={"-2": "Media", "-1": "Special", "0": "", "1": "Talk", "2": "User", "3": "User talk", "4": "Bioinformatics Core Wiki", "5": "Bioinformatics Core Wiki talk", "6": "File", "7": "File talk", "8": "MediaWiki", "9": "MediaWiki talk", "10": "Template", "11": "Template talk", "12": "Help", "13": "Help talk", "14": "Category", "15": "Category talk"},
wgNamespaceIds={"media": -2, "special": -1, "": 0, "talk": 1, "user": 2, "user_talk": 3, "bioinformatics_core_wiki": 4, "bioinformatics_core_wiki_talk": 5, "file": 6, "file_talk": 7, "mediawiki": 8, "mediawiki_talk": 9, "template": 10, "template_talk": 11, "help": 12, "help_talk": 13, "category": 14, "category_talk": 15, "image": 6, "image_talk": 7},
wgSiteName="Bioinformatics Core Wiki",
wgCategories=[],
wgRestrictionEdit=[],
wgRestrictionMove=[];
</script><script src="/skins/common/wikibits.js?270"></script>
<script src="/skins/common/ajax.js?270"></script>
<script src="/index.php?title=-&amp;action=raw&amp;gen=js&amp;useskin=monobook&amp;270"></script>

</head>
<body class="mediawiki ltr ns-0 ns-subject page-TrimBWAstyle_pl skin-monobook">
<div id="globalWrapper">
<div id="column-content"><div id="content" >
	<a id="top"></a>
	
	<h1 id="firstHeading" class="firstHeading">TrimBWAstyle.pl</h1>
	<div id="bodyContent">
		<h3 id="siteSub">From Bioinformatics Core Wiki</h3>
		<div id="contentSub"></div>
		<div id="jump-to-nav">Jump to: <a href="#column-one">navigation</a>, <a href="#searchInput">search</a></div>
		<!-- start content -->
<pre>

#!/usr/bin/perl

# AUTHOR: Joseph Fass
# LAST REVISED: January 2010
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2010 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;
use Getopt::Std;

my $usage = &quot;\nusage: cat original.fastq | perl $0 -q # &gt; trimmed.fastq\n\n&quot;.
            &quot;Trims fastq sequences a la Heng Li's bwa '-q' option.  &quot;.
            &quot;Sequences must be on one line each (no multi-line sequences).  &quot;.
            &quot;Sequences must be in Sanger fastq encoding ( phredScore = ord(Q)-33 )&quot;.
            &quot;Sequences trimmed down to zero length will have one base ('N') with quality 2 ('#').\n&quot;.
            &quot;-q #          targets # as quality level (default 20) ... NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)\n\n&quot;;
our($opt_q);
getopts('q:') or die $usage;
if (!defined($opt_q) or&nbsp;!($opt_q =~ m/^[0-9]+$/)) {$opt_q = 20;}

my $h1;  my $s;  my $h2;  my $q;
my $pos;  my $maxPos;  my $area;  my $maxArea;
while ($h1 = &lt;&gt;) {  # read first header
	$s = &lt;&gt;;  # read sequence
	chomp $s;
	$h2 = &lt;&gt;;  # read second header
	$q = &lt;&gt;;  # read quality scores
	chomp $q;
	$pos = length($q);
	$maxPos = $pos;
	$area = 0;
	$maxArea = 0;
	while ($pos&gt;0 and $area&gt;=0) {
		$area += $opt_q - (ord(substr($q,$pos-1,1))-33);
		if ($area &gt; $maxArea) {
			$maxArea = $area;
			$maxPos = $pos;
		}
		$pos--;
	}
	if ($pos==0) { $s = &quot;N\n&quot;;  $q = &quot;#\n&quot; }  # scanned whole read and didn't integrate to zero?  replace with &quot;empty&quot; read ...
	else {  # integrated to zero?  trim before position where area reached a maximum (~where string of qualities were still below 20 ...)
		$s = substr($s,0,$maxPos).&quot;\n&quot;;
		$q = substr($q,0,$maxPos).&quot;\n&quot;;
	}
	print $h1.$s.$h2.$q;
}

</pre>

<!-- 
NewPP limit report
Preprocessor node count: 4/1000000
Post-expand include size: 0/2097152 bytes
Template argument size: 0/2097152 bytes
Expensive parser function count: 0/100
-->

<!-- Saved in parser cache with key bioinf_mwiki:pcache:idhash:321-0!1!0!!en!2!edit=0 and timestamp 20110905174423 -->
<div class="printfooter">
Retrieved from "<a href="http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl">http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl</a>"</div>
		<div id='catlinks' class='catlinks catlinks-allhidden'></div>		<!-- end content -->
				<div class="visualClear"></div>
	</div>
</div></div>
<div id="column-one">
	<div id="p-cactions" class="portlet">
		<h5>Views</h5>
		<div class="pBody">
			<ul>
				 <li id="ca-nstab-main" class="selected"><a href="/index.php/TrimBWAstyle.pl" title="View the content page [c]" accesskey="c">Page</a></li>
				 <li id="ca-talk" class="new"><a href="/index.php?title=Talk:TrimBWAstyle.pl&amp;action=edit&amp;redlink=1" title="Discussion about the content page [t]" accesskey="t">Discussion</a></li>
				 <li id="ca-viewsource"><a href="/index.php?title=TrimBWAstyle.pl&amp;action=edit" title="This page is protected.&#10;You can view its source [e]" accesskey="e">View source</a></li>
				 <li id="ca-history"><a href="/index.php?title=TrimBWAstyle.pl&amp;action=history" title="Past revisions of this page [h]" accesskey="h">History</a></li>
			</ul>
		</div>
	</div>
	<div class="portlet" id="p-personal">
		<h5>Personal tools</h5>
		<div class="pBody">
			<ul>
				<li id="pt-login"><a href="/index.php?title=Special:UserLogin&amp;returnto=TrimBWAstyle.pl" title="You are encouraged to log in; however, it is not mandatory [o]" accesskey="o">Log in</a></li>
			</ul>
		</div>
	</div>
	<div class="portlet" id="p-logo">
		<a style="background-image: url(/skins/common/images/gcLogo.png);" href="/index.php/Main_Page" title="Visit the main page"></a>
	</div>
	<script type="text/javascript"> if (window.isMSIE55) fixalpha(); </script>
	<div class='generated-sidebar portlet' id='p-navigation'>
		<h5>Navigation</h5>
		<div class='pBody'>
			<ul>
				<li id="n-Home"><a href="/index.php/Main_Page">Home</a></li>
				<li id="n-Request-Service"><a href="/index.php/Request_Service">Request Service</a></li>
				<li id="n-Resources"><a href="/index.php/Computing_Resources">Resources</a></li>
				<li id="n-People-.2F-Contact-Info"><a href="/index.php/Contact">People / Contact Info</a></li>
				<li id="n-FAQ"><a href="/index.php/FAQ">FAQ</a></li>
			</ul>
		</div>
	</div>
	<div class='generated-sidebar portlet' id='p-Community'>
		<h5>Community</h5>
		<div class='pBody'>
			<ul>
				<li id="n-Bioinformatics-Mailing-List"><a href="/index.php/Bioinformatics_List">Bioinformatics Mailing List</a></li>
				<li id="n-BTF:-Bioinformatics-Technonoly-Forum"><a href="/index.php/BTF">BTF: Bioinformatics Technonoly Forum</a></li>
				<li id="n-UC-Davis-Bioinformatics-Tools"><a href="/index.php/UCD_BioTools">UC Davis Bioinformatics Tools</a></li>
			</ul>
		</div>
	</div>
	<div id="p-search" class="portlet">
		<h5><label for="searchInput">Search</label></h5>
		<div id="searchBody" class="pBody">
			<form action="/index.php" id="searchform">
				<input type='hidden' name="title" value="Special:Search"/>
				<input id="searchInput" title="Search Bioinformatics Core Wiki" accesskey="f" type="search" name="search" />
				<input type='submit' name="go" class="searchButton" id="searchGoButton"	value="Go" title="Go to a page with this exact name if exists" />&nbsp;
				<input type='submit' name="fulltext" class="searchButton" id="mw-searchButton" value="Search" title="Search the pages for this text" />
			</form>
		</div>
	</div>
	<div class="portlet" id="p-tb">
		<h5>Toolbox</h5>
		<div class="pBody">
			<ul>
				<li id="t-whatlinkshere"><a href="/index.php/Special:WhatLinksHere/TrimBWAstyle.pl" title="List of all wiki pages that link here [j]" accesskey="j">What links here</a></li>
				<li id="t-recentchangeslinked"><a href="/index.php/Special:RecentChangesLinked/TrimBWAstyle.pl" title="Recent changes in pages linked from this page [k]" accesskey="k">Related changes</a></li>
<li id="t-specialpages"><a href="/index.php/Special:SpecialPages" title="List of all special pages [q]" accesskey="q">Special pages</a></li>
				<li id="t-print"><a href="/index.php?title=TrimBWAstyle.pl&amp;printable=yes" rel="alternate" title="Printable version of this page [p]" accesskey="p">Printable version</a></li>				<li id="t-permalink"><a href="/index.php?title=TrimBWAstyle.pl&amp;oldid=2150" title="Permanent link to this revision of the page">Permanent link</a></li>			</ul>
		</div>
	</div>
</div><!-- end of the left (by default at least) column -->
<div class="visualClear"></div>
<div id="footer">
	<div id="f-poweredbyico"><a href="http://www.mediawiki.org/"><img src="/skins/common/images/poweredby_mediawiki_88x31.png" height="31" width="88" alt="Powered by MediaWiki" /></a></div>
	<ul id="f-list">
		<li id="lastmod"> This page was last modified on 19 January 2010, at 21:55.</li>
		<li id="viewcount">This page has been accessed 1,608 times.</li>
		<li id="privacy"><a href="/index.php/Bioinformatics_Core_Wiki:Privacy_policy" title="Bioinformatics Core Wiki:Privacy policy">Privacy policy</a></li>
		<li id="about"><a href="/index.php/Bioinformatics_Core_Wiki:About" title="Bioinformatics Core Wiki:About">About Bioinformatics Core Wiki</a></li>
		<li id="disclaimer"><a href="/index.php/Bioinformatics_Core_Wiki:General_disclaimer" title="Bioinformatics Core Wiki:General disclaimer">Disclaimers</a></li>
	</ul>
</div>
</div>

<script>if (window.runOnloadHook) runOnloadHook();</script>
<!-- Served in 0.324 secs. --></body></html>
