package WormBase::API::Object::Cds;

use Moose;

extends 'WormBase::API::Object';
with 'WormBase::API::Role::Object';
with 'WormBase::API::Role::Position';
with 'WormBase::API::Role::Sequence';

use Bio::Graphics::Browser2::Markup;

=pod

=head1 NAME

WormBase::API::Object::Cds

=head1 SYNPOSIS

Model for the Ace ?Cds class.

=head1 URL

http://wormbase.org/species/cds

=cut

has 'method' => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build_method',
    );

sub _build_method {
    my ($self) = @_;
    my $object = $self->object;
    my $class = $object->class;
    my $method = $object->Method;
    my $details = $method->Remark if $method;
    return {
        description => "the method used to describe the $class",
        data => ($method || $details) ? {
            method => $method && "$method",
            details => $details && "$details",
        } : undef
    };
}

has 'type' => (
    is  => 'ro',
    lazy_build => 1,
   );

sub _build_type {
    my ($self) = @_;
    my $s = $self->object;
    # figure out where this sequence comes from
    # should rearrange in order of probability
    my $type;
    if ($self ~~ '@Locus') {
        $type = 'WormBase CDS';
    }
    elsif ($s->Coding) {
        $type = 'predicted coding sequence';
    }
    elsif ($s->get('cDNA')) {
        ($type) = $s->get('cDNA');
    }
    elsif (eval{_is_merged($s)}) {
        $type = 'merged sequence entry';
    }
    else {
        $type = $s->Properties(1);
    }
    $type ||= 'unknown';
    return $type;
}


#######################################
#
# CLASS METHODS
#
#######################################

#######################################
#
# INSTANCE METHODS
#
#######################################

############################################################
#
# The Overview widget
#
############################################################

# name {}
# Supplied by Role

# taxonomy {}
# Supplied by Role

# description { }
# Supplied by Role

# sequence_type {}
# Supplied by Role

# identity {}
# Supplied by Role

# method {}
# Supplied by Role

# remarks {}
# Supplied by Role

# laboratory { }
# Supplied by Role


sub partial{
    my ($self) = @_;
    my $cds = $self->object;

    my $data = defined $cds->get('Start_not_found') || defined $cds->get('End_not_found') ?
        sprintf("%s%s%s not found",
            defined $cds->get('Start_not_found') ? "start" : "",
            defined $cds->get('Start_not_found') && defined $cds->get('End_not_found') ? " and " : "",
            defined $cds->get('End_not_found') ? "end" : ""
        ) : undef;

    print "start not found\n" if defined $cds->get('Start_not_found') ;
    print "end not found\n" if defined $cds->get('End_not_found');

    return {
        description => "Whether the start or end of the CDS is found",
        data => $data
    };
}

# gene_history { }
# This mehtod will return a data structure containing the
# historical record of the dead gene originally associated with this
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/cds/JC8.10a/gene_history


sub gene_history {
    my $self = shift;
    my $object = $self->object;

    my @historical_gene = map { $self->_pack_obj($_) } $object->Gene_history;
    return { description => 'Historical record of the genes originally associated with this CDS',
             data        => @historical_gene ? \@historical_gene : undef,
    };
}


############################################################
#
# The External Links widget
#
############################################################

# xrefs {}
# Supplied by Role

############################################################
#
# The Location Widget
#
############################################################

# genomic_position { }
# Supplied by Role

# sub tracks {}
# Supplied by Role

sub _build_tracks {
    my ($self) = @_;

    return {
        description => 'tracks to display in GBrowse',
        data => $self->_parsed_species =~ /elegans/
	    ? ($self->method->{data}{method} eq 'history') ?  [qw(HISTORICAL_GENES)] : [qw(GENES TRANSPOSONS TRANSPOSON_GENES EST_BEST PROTEIN_MOTIFS)] : undef,
    };
}

# genomic_image { }
# Supplied by Role

# note for AD:
# this one needs some reworking. it currently fetches the first segment
# in $self->segments, recomputes the start & stop, fetches more segments
# if the seq is a CDS or Transcript, and if more than 1 seg, selects the
# first one that matches the start and stop (or just the first one).
# throwing that segment back into genomic_position just wraps it up
sub _build_genomic_image {
    my ($self) = @_;
    my $seq = $self->object;
    return {
        description => 'The genomic location of the sequence to be displayed by GBrowse',
        data => undef
    } unless(defined $self->_segments && defined $self->_segments->[0] && $self->_segments->[0]->length< 100_0000);

    my $source = $self->_parsed_species;
    my $segment = $self->_segments->[0];

    my $ref   = $segment->ref;
    my $start = $segment->start;
    my $stop  = $segment->stop;

    # add another 10% to left and right
    $start = int($start - 0.1*($stop-$start));
    $stop  = int($stop  + 0.1*($stop-$start));
    my @segments;
    my $gene = eval { $seq->Gene;} || $seq;
    @segments = $self->gff->segment($gene);
    @segments = grep {$_->method eq 'wormbase_cds'} $self->gff->fetch_group(CDS => $seq)
	unless @segments;   # CB discontinuity
    # In cases where more than one segment is retrieved
    # (ie with EST or OST mappings)
    # choose that which matches the original segment.
    # This is slightly bizarre but expedient fix.
    my $new_segment;
    if (@segments > 1) {
        foreach (@segments) {
            if ($_->start == $start && $_->stop == $stop) {
                $new_segment = $_;
                last;
            }
        }
    }

    my ($position) = $self->_genomic_position([$new_segment || $segment || ()]);
    return {
        description => 'The genomic location of the sequence to be displayed by GBrowse',
        data        => $position,
    };
}

# genetic_position {}
# Supplied by Role



############################################################
#
# The Reagents Widget
#
############################################################

# orfeome_assays {}
# Supplied by Role

# microarray_assays {}
# Supplied by Role

# pcr_products {}
# Supplied by Role

# matching_cdnas {}
# Supplied by Role

# source_clone {}
# Supplied by Role

############################################################
#
# The Sequence Widget
#
############################################################

# print_blast {}
# Supplied by Role

# print_sequence {}
# Supplied by Role

# print_homologies {}
# Supplied by Role

# print_feature {}
# Supplied by Role

# strand {}
# Supplied by Role

# transcripts {}
# Supplied by Role

# predicted_units {}
# Supplied by Role

# predicted_exon_structure { }
# This method will return a data structure listing
# the exon structure contained within the sequence.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/cds/JC8.10a/predicted_exon_structure

sub predicted_exon_structure {
    my ($self) = @_;
    my $s = $self->object;

    my $index = 1;
    my @exons = map {
		my ($es,$ee) = $_->row;
		{
			no		=> $index++,
			start	=> "$es" || undef,
			end		=> "$ee" || undef,
			len 	=> "$es" && "$ee" ? $ee-$es+1 : undef
		};
	} $s->get('Source_Exons');

    return { description => 'predicted exon structure within the sequence',
             data        => @exons ? \@exons : undef };
}

############################################################
#
# PRIVATE METHODS
#
############################################################

sub _is_gap {
    return shift =~ /(\b|_)GAP(\b|_)/i;
}

sub _is_merged {
    return shift =~ /LINK|CHROMOSOME/i;
}

sub _build__segments {
    my ($self) = @_;
    my $object = $self->object;
    return [] unless $self->gff;
    # special case: return the union of 3' and 5' EST if possible
    if ($self->type =~ /EST/) {
        if ($object =~ /(.+)\.[35]$/) {
            my $base = $1;
            my ($seg_start) = $self->gff->segment("$base.3");
            my ($seg_stop)  = $self->gff->segment("$base.5");
            if ($seg_start && $seg_stop) {
                my $union = $seg_start->union($seg_stop);
                return [$union] if $union;
            }
        }
    }
    return [map {$_->absolute(1);$_} sort {$b->length<=>$a->length} $self->gff->segment($object)];
}




sub _print_unspliced {
    my ($self,$seq_obj,$unspliced,@features) = @_;
    my $name = $seq_obj->info . ' (' . $seq_obj->start . '-' . $seq_obj->stop . ')';

    my $length   = length $unspliced;
    if ($length > 0) {
        # mark up the feature locations

        my @markup;
        my $offset = $seq_obj->start;
        my $counter = 0;
        for my $feature (@features) {
            my $start    = $feature->start - $offset;
            my $length   = $feature->length;
            my $style = $feature->method eq 'CDS'  ? 'cds'.$counter++%2
            : $feature->method =~ /exon/ ? 'cds'.$counter++%2
            : $feature->method =~ 'UTR' ? 'utr' : '';
            push @markup,[$style,$start,$start+$length];
            push @markup,['uc',$start,$start+$length] unless $style eq 'utr';
        }
        push @markup,map {['space',10*$_]}   (1..length($unspliced)/10);
        push @markup,map {['newline',80*$_]} (1..length($unspliced)/80);
#       my $download = _to_fasta("$name|unspliced + UTR - $length bp",$unspliced);
        $self->markup_scheme->markup(\$unspliced,\@markup);
        return {
            #download => $download,
            header=>"unspliced + UTR",
            sequence=>$unspliced,
            length => $length,
            style=> 1,

        };
    }
}

# Fetch and markup the spliced DNA
# markup alternative exons
sub _print_spliced {
    my ($self, @features) = @_;
    my $spliced = join('',map {$_->dna} @features);
    my $splen   = length $spliced;
    my $last    = 0;
    my $counter = 0;
    my @markup  = ();
    my $prefasta = $spliced;
    for my $feature (@features) {
        my $length = $feature->stop - $feature->start + 1;
print "length!! $length\n";
        my $style  = $feature->method =~ /UTR/i ? 'utr' : 'cds' . $counter++ %2;
        my $end = $last + $length;
        push @markup,[$style,$last,$end];
        push @markup,['uc',$last,$end] if $feature->method =~ /exon/;
        $last += $length;
    }

    push @markup,map {['space',10*$_]}   (1..length($spliced)/10);
    push @markup,map {['newline',80*$_]} (1..length($spliced)/80);
    my $name = eval { $features[0]->refseq->name } ;
#   my $download=_to_fasta("$name|spliced + UTR - $splen bp",$spliced);
    $self->markup_scheme->markup(\$spliced,\@markup);

    return {                    # download => $download ,
        header=>"spliced + UTR",
        sequence=>$spliced,
        length=> $splen,
        style=> 1,
    } if $name;

}

sub _print_protein {
    my ($self,$features,$genetic_code) = @_;
#   my @markup;
    my $trimmed = join('',map {$_->dna} grep {$_->method eq 'coding_exon'} @$features);
    return unless $trimmed;     # Hack for mRNA
    my $peptide = Bio::Seq->new(-seq=>$trimmed)->translate->seq;
    my $change  = $peptide =~/\w+\*$/ ? 1 : 0;
    my $plen = length($peptide) - $change;

#   @markup = map {['space',10*$_]}      (1..length($peptide)/10);
#   push @markup,map {['newline',80*$_]} (1..length($peptide)/80);
    my $name = eval { $features->[0]->refseq->name };
#   my $download=_to_fasta("$name|conceptual translation - $plen aa",$peptide);
#   $markup->markup(\$peptide,\@markup);
    $peptide =~ s/^\s+//;

    return {                    # download => $download,
        header=>"conceptual translation",
        sequence=>$peptide,
        type => "aa",
        length => $plen,
    };
}

##use this or template to format sequence?

sub _to_fasta {
    my ($name,$dna) = @_;
    $dna ||= '';
    my @markup;
    for (my $i=0; $i < length $dna; $i += 10) {
        push (@markup,[$i,$i % 80 ? ' ':"\n"]);
    }
    _markup(\$dna,\@markup);
    $dna =~ s/^\s+//;
    $dna =~ s/\*$//;
    return  {   header=>"Genomic Sequence",
                content=>"&gt;$name\n$dna"
               };
}

# insert HTML tags into a string without disturbing order
sub _markup {
    my $string = shift;
    my $markups = shift;
    for my $m (sort {$b->[0]<=>$a->[0]} @$markups) { #insert later tags first so position remains correct
        my ($position,$markup) = @$m;
        next unless $position <= length $$string;
        substr($$string,$position,0) = $markup;
    }
}
# get coordinates of parent for exons etc
sub _get_parent_coords {
    my ($self,$s) = @_;
    my ($parent) = $self->sequence;
    return unless $parent;
    #  my $subseq = $parent->get('Subsequence');  # prevent automatic dereferencing

    # Escape the sequence name for fetching
    $s =~ s/\./\\./g;
    # We may be dealing with transcripts, too.
    my $se;
    foreach my $tag (qw/CDS_child Transcript/) {
        my $subseq = $parent->get($tag); # prevent automatic dereferencing
        if ($subseq) {
            $se = $subseq->at($s);
            if ($se) {
                my ($start,$stop) = $se->right->row;
                my $orientation = $start <=> $stop;
                return ($start,-$orientation,$parent);
            }
        }
    }
    return;
}

__PACKAGE__->meta->make_immutable;

1;
