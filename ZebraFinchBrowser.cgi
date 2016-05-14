#!/usr/bin/perl

use strict;
use warnings;
use CGI qw/:standard/;
#use CGI::Carp qw/ fatalsToBrowser warningsToBrowser/;
use Bio::Graphics;
use Bio::SeqFeature::Generic;

use lib '/home/john.samuel/src/ensembl/modules';
print "Content-type: text/html\n\n";

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# if information has been entered into the form and submitted
if (param()) {
    # take entered values and store to variables
    my $chromosome = param('chromosome');
    my $start = param('start');
    my $end = param('end');
    my $species = "Taeniopygia guttata";
    my @errors;
    
    print begin_html("Report for Zebra Finch $chromosome: $start-$end", "Zebra Finch Genome Browser");
    
    # validation from form  
    if ($start !~ /^\d+$/g ) {
        push @errors, "Start position must only contain digits";
    }
    if ($end !~ /^\d+$/g ) {
        push @errors, "End position must only contain digits";
    }
    if ($start < 1) {
        push @errors, "Start position cannot be less than 1";
    }
    
    # creates object of zebra finch species
    my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
    # creates a slice of the user selected chromosome, used to validate end position later
    my $chromo_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome);
    # stores the end position of the chromosome selected by user (for validation)
    my $max_end = $chromo_slice->end();
    
    #more validation using information from chromosome object
    if ($end > $max_end) {
        push @errors, "The end position rests outside the selected chromosome";
    }
    if ($start > $end) {
        push @errors, "The start position must be less than the end position";
    }
    if ($end - $start < 1e3) {
        push @errors, "The distance between start and end positions must be greater than 1,000bp";
    }
    if ($end - $start > 10e7) {
        push @errors, "The distance between start and end positions must be less than 10,000,000bp";
    }
    # if any of the validation steps fail, errors pushed into this array
    # this prints all the acquired errors, if any
    if (@errors) {
        #reprints instructions with examples
        print instructions();
        foreach my $error (@errors) {
            print "<font color='red'>$error</font><br>";
        }
        # reprints the form, with previous (error prone) entries saved
        print form($chromosome, $start, $end);
        # exits program, otherwise creates errors from ensembl 
        exit;
        print end_html();
    }
    #prints the report information if valid information is supplied by the user
    else {
        print "Report for Zebra Finch $chromosome: $start-$end.<br>";
    }
    
    my $line_length = $end - $start;
    # creates panel for the diagram to become 'attached' to
    my $panel = Bio::Graphics::Panel->new(-length => $line_length, -width  => 800, -pad_top=> 50, -pad_left=>100, -pad_right=>100,
		-start=>$start,-end=>$end);
    # creates the line/scale based on the start and end positions supplied by user
    my $line = Bio::SeqFeature::Generic->new(-start => $start, -end=> $end);
    # adds the line/scale to the panel
    $panel->add_track($line, -glyph => 'arrow',-tick => 2,-fgcolor => 'black',-double  => 1);
    
    # creates an object of the entire chromosome selected
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, $start, $end);
    # fetches all the genes located within the selected chromosome
    my @genes = @{ $slice->get_all_Genes() };
    
    # prints the table header with information about each gene
    print <<T_HEADER;
    <table border ="1">
    <tr><th>Gene ID</th><th>Start Position</th><th>End Position</th><th>Strand</th><th>Length</th><th>Description</th><th>External Name</th><th>Gene Type</th><th>Status</th><th>Number of Transcripts</td></tr>
T_HEADER

    #perform these tasks for each of the genes located within this region, if there are genes in this given region
    if (@genes) {
        foreach my $gene (@genes) {
            my $gene_id = $gene->stable_id();
            my $start = $gene->seq_region_start();
            my $end = $gene->seq_region_end();
            my $strand = $gene->strand();
            my $length = $gene->length();
            
            # pre tags maintain the html whitespace so border lines are displayed even if description and external name are null 
            my $description = "<pre> </pre>";
            
            if ($gene->description()) {
                $description = $gene->description();
            }
             
            my $external_name = "<pre> </pre>";
            
            if($gene->external_name()) {
                $external_name = $gene->external_name();
            }
            my $gene_type = $gene->biotype();
            my $status = $gene->status();
            
            #grabs all the transcripts within each given gene, stores them to an array
            my @transcript = @{ $gene->get_all_Transcripts };
            # counts the number of transcripts within each gene
            my $transcript_num = scalar(@transcript);
            my $color;
            
            #sets the color of the diagram biotype font depending on what gene type it is.
            if ($gene_type eq "protein_coding") {
                $color = "red";
            }
            else {
                $color = "black";
            }
            
            # adds each gene within the designated region to the diagram
            my $track = $panel->add_track(-glyph => 'transcript2', -stranded =>1, -label => 1, -fontcolor => $color, 
                -bgcolor => 'green', -description=>$gene_type);
        
            my $diagram_info = "$gene_id: ($start-$end)";
            
            # creates an object that has the information about the gene
            my $display = Bio::SeqFeature::Generic->new(-display_name => $diagram_info, -start => $start, -end => $end, -strand=>$strand);
            #adds the information for each gene to the diagram
            $track->add_feature($display);
            
            # subsitutes 1/-1 to +/- for the table
            if ($strand == 1) {
                $strand =~ s/1/+/;
            }
            if ($strand == -1) {
                $strand =~ s/-1/-/;
            }
            
            # prints the information of every gene within the given region into the table
            print "<div class ='test'><tr><td><a href='http://uswest.ensembl.org/Gene/Summary?db=core;g=$gene_id' target='_blank'>$gene_id</td><td>$start</td><td>$end</td><td>$strand</td><td>$length</td><td>$description</td>
            <td>$external_name</td><td>$gene_type</td><td>$status</td><td>$transcript_num</td></tr></div>";
        }
    print "</table>";
    
    # creates a png file
    open(my $fh, ">", "browser.png") or die "cannot open > browser.png: $!";
    # prints the image to the file
    print $fh $panel->png;
    close $fh;
    # add image/diagram to the webpage
    print "<center><img src='browser.png'/></center>";
    }
    
    #if the selected region does not contain any genes, print this inside the table
    else {
        print "<tr><td colspan ='10' align='middle'>Zero genes found</tr>";
        print "</table>";
    }
    # link to go back and enter another search
    print "<center><a href='http://zenit.senecac.on.ca/~bif724_161a18/Assignment3/tash1.cgi'>New Search</a></center>";
    print bottom_html();
}

#if the form has not been submitted, print the instructions and form (contaisn default values)
else {
    #form not submitted so print form for user
    #print "Content-type: text/html\n\n";
    print begin_html("Select Region", "Zebra Finch Genome Browser");
    print instructions();
    print form("1A", "38697201", "39070179");
    print bottom_html();
}


# subroutines
#-------------------------------
sub begin_html {
    my $tabHeading = shift;
    my $h1 =shift;
    return<<begin_html;
<!DOCTYPE html>
<html>
    <head>
        <title>$tabHeading</title>
    </head>
    <body>
    <h1 align='center'>$h1</h1>
begin_html
}

sub bottom_html {
    return<<bottom_html;
    </body>
</html>
bottom_html
}

sub form {
    #my $form = "";
    #my $chromosome = param("chromosome");
    #my $start = param('start');
    #my $end = param('end');
    my $form = "";
    my $chromosome = shift;
    my $start = shift;
    my $end = shift;
    
    
    #print "<h1 align='center'>Zebra Finch Genome Browser</h1>";
  
    $form .= <<FORM;
        <form action="$0" method="get">
            <div class="table">
                <table class="add">
                <tr><td>Choose chromosome</td><td>
                    <select name ="chromosome"/>
FORM
    my @chromosomes = qw/ 1 1A 1B 2 3 4 4A 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 LGE22 LG2 LG5 Z MT/;
    foreach (@chromosomes) {
        my $selection = "";
        $selection = " selected='selected'" if ($chromosome eq "$_");
        $form .= "<option $selection>$_</option>";
    }
    $form .=<<FORM;
                    </select>
                </tr></td>
                <tr><td>Start position</td><td><input type ="text" class="textbox" name ="start" value ="$start"</td></tr>
                <tr><td>End position</td><td><input type ="text" class="textbox" name ="end" value ="$end"</td></tr>
                </td></tr><tr><td></td><td colspan = "2" align = "right"><input type ="submit"></td></tr> 
                </table>
            </div>
        </form>
FORM
    return $form;
}
sub instructions  {
    return<<INSTRUCT;
        <center><p>The zebra finch (<i>Taeniopygia guttata</i>), native to Australia is known for its distinctive calls, "beep", "meep", and even "a-ha"! Use this
        customization of the ensembl genome browser to view regions of its genome! Just enter the desired chromosome and the start and end positions. Some example regions to explore:</p></center>
        
        <table border="1">
            <tr><th>Chromosome</th><th>Start Position</th><th>End Position</th><th>Number of Genes</th><tr>
            <tr><td>1A</td><td>38697201</td><td>39070179</td><td>4</td></tr>
            <tr><td>1</td><td>1</td><td>11000</td><td>0</td></tr>
            <tr><td>LGE22</td><td>1</td><td>25000</td><td>1</td></tr>
        </table>
        
            
        
INSTRUCT

}
