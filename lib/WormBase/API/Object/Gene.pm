
package WormBase::API::Object::Gene;

use Moose;
use File::Spec::Functions qw(catfile catdir);
use namespace::autoclean -except => 'meta';
use File::Temp;
use JSON;
use HTTP::Tiny;
use Data::Dumper;

extends 'WormBase::API::Object';
with    'WormBase::API::Role::Object';
with    'WormBase::API::Role::Position';
with    'WormBase::API::Role::Interaction';
with    'WormBase::API::Role::Variation';
with    'WormBase::API::Role::Expression';
with    'WormBase::API::Role::Feature';

=pod

=head1 NAME

WormBase::API::Object::Gene

=head1 SYNPOSIS

Model for the Ace ?Gene class.

=head1 URL

http://wormbase.org/species/*/gene

=head1 METHODS/URIs

=cut

has '_all_proteins' => (
    is  => 'ro',
    lazy => 1,
    default => sub {
        return [
            map { $_->Corresponding_protein(-fill => 1) }
                shift->object->Corresponding_CDS
        ];
    }
);

has 'sequences' => (
    is  => 'ro',
    lazy => 1,
    builder => '_build_sequences',
);

sub _build_sequences {
    my $self = shift;
    my $gene = $self->object;
    my %seen;
    my @seqs = grep { !$seen{$_}++} $gene->Corresponding_transcript;

    for my $cds ($gene->Corresponding_CDS) {
        next if defined $seen{$cds};
        my @transcripts = grep {!$seen{$cds}++} $cds->Corresponding_transcript;

        push (@seqs, @transcripts ? @transcripts : $cds);
    }
    return \@seqs if @seqs;
    return [$gene->Corresponding_Pseudogene];
}

has 'tracks' => (
    is      => 'ro',
    lazy    => 1,
    default => sub {
        my $self = shift;
        return {
            description => 'tracks displayed in GBrowse',
            data        => $self->object->Corresponding_transposon ? [qw/TRANSPOSONS TRANSPOSON_GENES/] : [qw/GENES VARIATIONS_CLASSICAL_ALLELES CLONES/],
        };
    }
);

has '_alleles' => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build__alleles',
);

sub _build__alleles {
    my ($self) = @_;
    my $object = $self->object;

    my $count = $self->_get_count($object, 'Allele');
    my @all = $object->Allele;
    my @alleles;
    my @polymorphisms;

    foreach my $allele (@all) {
              (grep {/Natural_variant|RFLP/} $allele->Variation_type) ?
                    push(@polymorphisms, $self->_process_variation($allele)) :
                    push(@alleles, $self->_process_variation($allele));
    }

    return {
        alleles        => @alleles ? \@alleles : undef,
        polymorphisms  => @polymorphisms ? \@polymorphisms : undef,
    };

}

# sub phenotype_by_interaction {}

# Delegated to REST.

#######################################
#
# The Overview Widget
#   template: classes/gene/overview.tt2
#
#######################################

# also_refers_to { }
# This method will return a data structure containing
# other names that have also been used to refer to the
# gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/also_refers_to

# Delegate to REST

# named_by { }
# This method will return a data structure containing
# the WB person who named the gene
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/named_by

# Delegate to REST

# classification { }
# This method will return a data structure containing
# the general classification of the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/classification

# Delegated to REST

# cloned_by { }
# This method will return a data structure containing
# the person or laboratory who cloned the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/cloned_by

# Delegated to REST.

# parent_sequence { }
# This method will return a data structure containing
# the parent sequence of the gene
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/parent_sequence

# Delegated to REST.

# clone { }
# This method will return a data structure containing
# the parent clone of the gene
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/clone

# Delegated to REST.

# concise_desciption { }
# This method will return a data structure containing
# the prose concise description of the gene, if one exists.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/concise_description

# Delegated to REST.

# gene_class { }
# This method will return a data structure containing
# the gene class packed tag of the gene, if one exists.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/gene_class

# Delegated to REST.

# operon { }
# This method will return a data structure containing
# the operon packed tag of the gene, if one exists.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/operon

# Delegated to REST

# transposon { }
# This method will return a data structure containing
# the transposon packed tag of the gene, if one exists.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/operon

# Delegated to REST

# legacy_information { }
# This method will return a data structure containing
# legacy information from the original Cold Spring Harbor
# C. elegans I & II texts.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/legacy_information

# Delegated to REST

# locus_name { }
# This method will return a data structure containing
# the name of the genetic locus.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/locus_name

# Delegated to REST

# name {}
# Supplied by Role

# other_names {}
# Supplied by Role

# sequence_name { }
# This method will return a data structure containing
# the primary sequence name of the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/sequence_name

# Delegated to REST

# status {}
# Supplied by Role

# structured_description { }
# This method will return a data structure containing
# various structured descriptions of gene's function.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/structured_description

# Delegated to REST.

# human_disease_relevance { }
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/human_disease_relevance

# Delegated to REST.

# taxonomy {}
# Supplied by Role

# version { }
# This method will return a data structure containing
# the current WormBase version of the gene.
# curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/version

# Delegated to REST

# merged_into {}

# Delegated to REST.

#######################################
#
# The External Links widget
#   template: shared/widgets/xrefs.tt2
#
#######################################

# xrefs {}
# Supplied by Role

#######################################
#
# The Genetics Widget
#   template: classes/gene/genetics.tt2
#
#######################################

# alleles { }
# This method will return a complex data structure
# containing alleles of the gene (but not including
# polymorphisms or other natural variations.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/alleles

sub alleles {
    my ($self, $dummy, $data) = @_;
    return $data;
}

# polymorphisms { }
# This method will return a complex data structure
# containing polymorphisms and natural variations
# but not alleles.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/polymorphisms

sub polymorphisms {
    my ($self, $dummy, $data) = @_;
    return $data;
}

# reference_allele { }
# This method will return a complex data structure
# containing the reference allele of the gene, if
# one exists.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/reference_allele

# Delegated to REST.

# strains { }
# This method will return a complex data structure
# containing strains carrying the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/strains

sub strains {
    # Delegate to REST.  However, we need a method here to override strains in API::Object.
    my ($self, $dummy, $data) = @_;
    return $data;
}

# rearrangements { }
# This method will return a data structure
# containing rearrangements affecting the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/rearrangements

sub rearrangements {
    my $self    = shift;
    my $object  = $self->object;
    my @positive = map { $self->_pack_obj($_) } $object->Inside_rearr;
    my @negative = map { $self->_pack_obj($_) } $object->Outside_rearr;

    return { description => 'rearrangements involving this gene',
             data        => (@positive || @negative) ? { positive => \@positive,
                             negative => \@negative
            } : undef
    };
}


#######################################
#
# The Gene Ontology widget
#   template: classes/gene/gene_ontology.tt2
#
#######################################

# gene ontology { }
# This method will return a data structure containing
# curated and electronically assigned gene ontology
# associations.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/gene_ontology

has '_gene_ontology' => (
    is  => 'ro',
    lazy => 1,
    builder => '_build_gene_ontology',
);
use Data::Dumper;
sub _build_gene_ontology {
    my $self   = shift;
    my $object = $self->object;

    my @data;
    foreach my $anno ($object->GO_annotation) {
        my $go_term = $anno->GO_term;
        my $go_type = $go_term->Type;
        my $go_code = $anno->GO_code;
        my $relation = $anno->Annotation_relation;

        my @entities = map {
            $self->_pack_list([$_->col()]);
        } $anno->Annotation_made_with;

        my %extensions = map {
            my ($ext_type, $ext_name, $ext_value) = $_->row();
           "$ext_name" => $self->_pack_obj($ext_value)
        } $anno->Annotation_extension;

        my $ev_names = ['Reference', 'Contributed_by', 'Date_last_updated'];
        my $evidence = $self->_get_evidence($anno->fetch(), $ev_names);

        my @term_details = () ; #('' . $go_term->Term);
        push @term_details, { evidence => \%extensions } if %extensions;

        my $anno_data = {
            term_id => $self->_pack_obj($go_term, "$go_term"),
            term_description => $self->_pack_obj($go_term), #'' . $go_term->Term,
            anno_id => "$anno",
            term => @term_details ? \@term_details : undef,
            evidence_code => $evidence ? { evidence => $evidence, text => "$go_code" } : "$go_code",
            go_type => "$go_type",
            with => @entities ? \@entities : undef,
            extensions => %extensions ? \%extensions : undef,
        };

        push @data, $anno_data;
    }

    return \@data;
}

sub gene_ontology {
    my ($self)   = @_;

    my %data_by_type = ();
    my @data = @{ $self->_gene_ontology };
    foreach my $anno_data (@data){
        my $type = $anno_data->{go_type};
        $data_by_type{$type} ||= ();
        push @{$data_by_type{$type}}, $anno_data;
    }

    return {
        description => 'gene ontology assocations',
        data        => %data_by_type ? \%data_by_type : undef,
    };
}

sub gene_ontology_summary {
    my ($self)   = @_;

    my @data = @{ $self->_gene_ontology };

    sub _get_type {
        my ($item) = @_;
        return $item->{go_type};
    }

    sub _get_go_term {
        my ($item) = @_;
        return $item->{term_id}->{label};
    }

    my $data_by_type = $self->_group_and_combine(\@data, \&_get_type);
    my %result_by_type = ();

    foreach my $go_type (keys %$data_by_type) {
        my $data4type = $data_by_type->{$go_type};
        my $result4type = $self->_group_and_combine($data4type, \&_get_go_term, \&_summarize_go_term);
        $result_by_type{$go_type} = [values %$result4type];
    }
    print Dumper \%result_by_type;
    return {
        description => 'gene ontology assocations',
        data        => %result_by_type ? \%result_by_type : undef,
    };
}

sub _group_and_combine {
    my ($self, $all_items, $group_fun, $summarize_fun) = @_;
    my %grouped = ();

    # to group
    foreach my $item (@$all_items) {
        my $group_key = $group_fun->($item);
        $grouped{$group_key} ||= ();
        push  @{ $grouped{$group_key} }, $item;
    }

    # to summarize
    my %summarized = ();
    if ($summarize_fun){
        foreach my $group_key (keys %grouped) {
            my @group_items = @{ $grouped{$group_key} };
            my $group_result = $summarize_fun->(\@group_items);
            $summarized{$group_key} = $group_result;
        }
    }
    return %summarized ? \%summarized : \%grouped;
}
sub _summarize_go_term {
    my ($anno_data_all) = @_;

    my @exts_all = ();
    foreach my $anno_data (@$anno_data_all){
        #extensions within a single annotation
        my $exts = $anno_data->{extensions};
        push @exts_all, { evidence => $exts } if $exts;

        # my $exts = $anno_data->{extensions} || {};
        # my @pairs = map {[$_,$exts->{$_}]} (keys %$exts);
        # push @exts_all, \@pairs if @pairs;
    }
    return {
        term_id => $anno_data_all->[0]->{term_id},
        term_description => $anno_data_all->[0]->{term_description},
        extensions => @exts_all ? \@exts_all : undef,
    };
}

#######################################
#
# The History Widget
#    template: shared/widgets/history.tt2
#
#######################################

# history { }
# This method returns a data structure containing the
# curatorial history of the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/history

# Delegated to REST.

# Subroutine for the "Historical Annotations" table

# Delegated to REST.

#######################################
#
# The Homology Widget
#   template: classes/gene/homology.tt2
#
#######################################

# best_blastp_matches {}
# Supplied by Role

# nematode_orthologs { }
# This method returns a data structure containing the
# orthologs of this gene to other nematodes housed
# at WormBase.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/nematode_orthologs


# human_orthologs { }
# This method returns a data structure containing the
# human orthologs of this gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/human_orthologs

has '_other_orthologs' => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build__other_orthologs',
);

sub _build__other_orthologs {
    my ($self) = @_;
    return $self->_parse_homologs(
        [ $self->object->Ortholog_other ],
        sub {
            return $_[0]->right ? [map { $self->_pack_obj($_) } $_->right->col] : undef;
        }
    );
}

# I sure do wish we had some descriptions for human genes.


# other_orthologs { }
# This method returns a data structure containing the
# orthologs of this gene to species outside of the core
# nematodes housed at WormBase. See also nematode_orthologs()
# and human_orthologs();
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/other_orthologs


# paralogs { }
# This method returns a data structure containing the
# paralogs of this gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/paralogs

# Private helper method to standardize structure of homologs.
sub _parse_homologs {
    my ($self, $homologs, $method_sub) = @_;

    my @parsed;
    foreach (@$homologs) {
        my $packed_homolog = $self->_pack_obj($_);
        my $species = $packed_homolog->{taxonomy};
        my ($g, $spec) = split /_/, $species;
        push @parsed, {
            ortholog => $packed_homolog,
            method   => $method_sub->($_),
            species  => {
                genus   => ucfirst $g,
                species => $spec,
            },
        };
    }

    return \@parsed;
}

# human_diseases { }
# This method returns a data structure containing disease
# processes that human orthologs of this gene are thought
# to participate in.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/human_diseases

# Delegated to REST.

# protein_domains { }
# This method returns a data structure containing the
# protein domains contained in this gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/protein_domains

# Delegate to REST.

# treefam { }
# This method returns a data structure containing the
# link outs to the Treefam resource.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/treefam

sub treefam {
    my $self   = shift;
    my $object = $self->object;

    my %data;
    foreach (@{$self->_all_proteins}) {
        my $treefam = $self->_fetch_protein_ids($_,'treefam');
        # Ignore proteins that lack a Treefam ID
        next unless $treefam;
        $data{"$treefam"} = "";
    }
    my @data = keys %data;
    return { description => 'data and IDs related to rendering Treefam trees',
             data        => @data ? \@data : undef,
    };
}


#######################################
#
# The Location Widget
#
#######################################

# genomic_position { }
# Supplied by Role

sub _build_genomic_position {
    my ($self) = @_;
    my @pos = $self->_genomic_position([ $self->_longest_segment || () ]);
    return {
        description => 'The genomic location of the sequence',
        data        => @pos ? \@pos : undef,
    };
}

# genetic_position { }
# Supplied by Role

# sub genomic_image { }
# Supplied by Role


#######################################
#
# The Mapping Data Widget
#   template: classes/gene/mapping_data.tt2
#
#######################################

# These are currently kind-of inefficient because they're fetching the whole data
# for each field.  Can we add caching?

# two_pt_data
# this method returns mapping data associated with this gene
# adapted from old site
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/two_pt_data


# Delegated to REST.

# multi_pt_data
# this method returns mapping data associated with this gene
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/multi_pt_data

# Delegated to RESt.

# pos_neg_data
# this method returns mapping data associated with this gene
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/pos_neg_data

# Delegated to REST.

#######################################
#
# The Phenotype Widget
#
#######################################

# phenotype { }
# returns the phenotype(s) associated with the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/phenotype

# sub phenotype {
#     my ($self, %args) = @_;
#     my $data = $args{_data};
#     return $data;
# }



sub drives_overexpression {
    my ($self) = @_;
    my $object = $self->object;

    my %phenotypes;
    foreach my $construct ($object->Construct_product) {

        foreach my $transgene ($construct->Transgene_construct){

            my $summary = $transgene->Summary;

            # Retain in case we also want to add not_observed...
            foreach my $obs ('Phenotype'){
                foreach my $phene ($transgene->$obs){

                    # Only include those transgenes where the Caused_by in #Phenotype_info
                    # is the current gene.
                    my ($caused_by) = $phene->at('Caused_by');
                    next unless $caused_by eq $object;

                    $phenotypes{$obs}{$phene}{object} //= $self->_pack_obj($phene);
                    my $evidence = $self->_get_evidence($phene);
                    # $evidence->{Summary}   = "$summary" if $summary;
                    $evidence->{Transgene} = $self->_pack_obj($transgene);

                    if ($evidence && %$evidence){
                        my $transgene_label = $transgene->Public_name;
                        my $ev = {
                            text  => [$self->_pack_obj($transgene, $transgene_label),
                                      "<em>$summary</em>",
                                      $evidence->{Remark}],
                            evidence => $evidence
                        };
                        push @{$phenotypes{$obs}{$phene}{evidence}}, $ev;
                    }

                }
            }
        }
    }
    return { data        => (defined $phenotypes{Phenotype}) ? \%phenotypes : undef ,
             description => 'phenotypes due to overexpression under the promoter of this gene', };

}






#######################################
#
# The Reagents Widget
#
#######################################

# All done via REST.

#######################################
#
# The Regulation Widget
#   template: classes/gene/regulation.tt2
#
#######################################

# regulation_on_expression_level { }
# This method returns a data structure containing the
# a data table describing the regulation on expression
# level.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene000066763/regulation_on_expression_level

sub regulation_on_expression_level {
    my $self   = shift;
    my $object = $self->object;
    my $datapack = {
        description => 'Regulation on expression level',
        data        => undef,
    };
    return $datapack unless ($object->Gene_regulation);

    my @stash;

    # Explore the relationship in both directions.
    foreach my $tag (qw/Trans_regulator Trans_target/) {
        my $join = ($tag eq 'Trans_regulator') ? 'regulated by' : 'regulates';
        if (my @gene_reg = $object->$tag(-filled=>1)) {
            foreach my $gene_reg (@gene_reg) {
                my ($string,$target);
                if ($tag eq 'Trans_regulator') {
                    $target = $gene_reg->Trans_regulated_gene(-filled=>1)
                    || $gene_reg->Trans_regulated_seq(-filled=>1)
                    || $gene_reg->Other_regulated(-filled=>1);
                } else {
                    $target = $gene_reg->Trans_regulator_gene(-filled=>1)
                    || $gene_reg->Trans_regulator_seq(-filled=>1)
                    || $gene_reg->Other_regulator(-filled=>1);
                }
                # What is the nature of the regulation?
                # If Positive_regulate and Negative_regulate are present
                # in the same gene object, then it means the localization is changed.  Go figure.
                if ($gene_reg->Positive_regulate && $gene_reg->Negative_regulate) {
                    $string .= ($tag eq 'Trans_regulator')
                    ? 'Changes localization of '
                    : 'Localization changed by ';
                } elsif ($gene_reg->Result
                         and $gene_reg->Result eq 'Does_not_regulate') {
                    $string .= ($tag eq 'Trans_regulator')
                    ? 'Does not regulate '
                    : 'Not regulated by ';
                } elsif ($gene_reg->Positive_regulate) {
                    $string .= ($tag eq 'Trans_regulator')
                    ? 'Positively regulates '
                    : 'Positively regulated by ';
                } elsif ($gene_reg->Negative_regulate) {
                    $string .= ($tag eq 'Trans_regulator')
                    ? 'Negatively regulates '
                    : 'Negatively regulated by ';
                }

                # _pack_obj may already take care of this:
                my $common_name = $self->_public_name($target);
                push @stash, {
                    string          => $string,
                    target          => $self->_pack_obj($target, $common_name),
                    gene_regulation => $self->_pack_obj($gene_reg)
                };
            }
        }
    }

    $datapack->{data} = \@stash if @stash;
    return $datapack;
}

#######################################
#
# The References Widget
#
#######################################

# references {}
# Supplied by Role

#######################################
#
# The Sequences Widget
#
#######################################

# gene_models { }
# This method will return an extensive data structure containing
# gene models for the gene.
# eg: curl -H content-type:application/json http://api.wormbase.org/rest/field/gene/WBGene00006763/gene_models

# Delegated to REST.

# TH: Retired 2011.08.17; safe to delete or transmogrify to some other function.
# should we return entire sequence obj or just linking/description info? -AC
sub other_sequences {
    my $self   = shift;

    my @data = map {
        my $title = $_->Title;
        {
            sequence => $self->_pack_obj($_),
            description => $title && "$title",
        }
    } $self->object->Other_sequence;

    return {
        description => 'Other sequences associated with gene',
        data        => @data ? \@data : undef,
    };
}

#######################################
#
# The "Sequence features" Widget
#
#######################################

# features {}
# Supplied by Role

# Display a gbrowse image specific for the
# Sequence features widget.
sub feature_image {
    my ($self) = @_;

    my $segment = $self->_longest_segment;
    return unless $segment;

    # Create a NEW segment from this with expanded coordinates.
    my $dbh = $self->gff_dsn();# || return \@segments;
    my $start = $segment->start - 2000;
    my $stop  = $segment->stop  + 2000;
    my ($expanded_segment) = $dbh->segment($segment->seq_id,$start,$stop);

    my $position = $self->_seg2posURLpart($expanded_segment);
    $position->{tracks} = [qw/GENES RNASEQ_ASYMMETRIES RNASEQ RNASEQ_SPLICE POLYSOMES MICRO_ORF DNASEI_HYPERSENSITIVE_SITE REGULATORY_REGIONS PROMOTER_REGIONS HISTONE_BINDING_SITES TRANSCRIPTION_FACTOR_BINDING_REGION TRANSCRIPTION_FACTOR_BINDING_SITE BINDING_SITES_PREDICTED BINDING_SITES_CURATED BINDING_REGIONS/];

    return {
        description => 'The genomic location of the sequence to be displayed by GBrowse',
        data        => $position,
    };
}



#########################################
#
#   INTERNAL METHODS
#
#########################################

# This is for GO processing
# TH: I don't understand the significance of the nomenclature.
# Oh wait, I see, it's used to force an order in the view.
# This should probably be an attribute or view configuration.
sub _go_method_detail {
    my ($self,$method,$detail) = @_;
    return 'a_Curated' if $method =~ m/Paper/;
    return 'z_No Method' unless $detail;
    return 'b_Phenotype to GO Mapping' if ($detail =~ m/phenotype/i);
    return 'c_Interpro to GO Mapping' if ($detail =~ m/interpro/i);
    return 'd_TMHMM to GO Mapping' if ($detail =~ m/tmhmm/i);
    return 'z_No Method';
}

# Fetch unique transcripts (Transcripts or Pseudogenes) for the gene
sub _fetch_transcripts { # pending deletion
    my $self = shift;
    my $object = $self->object;
    my %seen;
    my @seqs = grep { !$seen{$_}++} $object->Corresponding_transcript;
    my @cds  = $object->Corresponding_CDS;
    foreach (@cds) {
        next if defined $seen{$_};
        my @transcripts = grep {!$seen{$_}++} $_->Corresponding_transcript;
        push (@seqs,(@transcripts) ? @transcripts : $_);
    }
    @seqs = $object->Corresponding_Pseudogene unless @seqs;
    return \@seqs;
}

sub _build__segments {
    my ($self) = @_;
    my $sequences = $self->sequences;
    my @segments;

    my $dbh = $self->gff_dsn();# || return \@segments;

    my $object = $self->object;

    if (@segments = $dbh->segment($object)
        ||  map {$dbh->segment($_)} @$sequences
        ||  map { $dbh->segment($_) } $object->Corresponding_Pseudogene # Pseudogenes (B0399.t10)
        ||  map { $dbh->segment( $_) } $object->Corresponding_Transcript # RNA transcripts (lin-4, sup-5)
    ) {
        return defined $segments[0] ? \@segments : undef;
    }

    return;
}

sub _build__gene {
    my ($self) = @_;
    my $object = $self->object;

    return $object;
}

# TODO: Logically this might reside in Model::GFF although I don't know if it is used elsewhere
# Find the longest GFF segment
sub _longest_segment {
    my ($self) = @_;
    # Uncloned genes will NOT have segments associated with them.
    my ($longest)
        = sort { $b->stop - $b->start <=> $a->stop - $a->start}
    @{$self->_segments} if $self->_segments;

    return $longest;
}

sub _select_protein_description { # pending deletion
    my ($self,$seq,$protein) = @_;
    my %labels = (
        Pseudogene => 'Pseudogene; not attached to protein',
        history     => 'historical prediction',
        RNA         => 'non-coding RNA transcript',
        Transcript  => 'non-coding RNA transcript',
    );
    my $error = $labels{eval{$seq->Method}};
    $error ||= eval { ($seq->Remark =~ /dead/i) ? 'dead/retired gene' : ''};
    my $msg = $protein ? $protein : $error;
    return $msg;
}


# I need to retain this in order to link to Treefam.
sub _fetch_protein_ids {
    my ($self,$s,$tag) = @_;
    my @dbs = $s->Database;
    foreach (@dbs) {
        return $_->right(2) if (/$tag/i);
    }
    return;
}

# TODO: This could logically be moved into a template
sub _other_notes { # pending deletion
    my ($self,$object) = @_;

    my @notes;
    if ($object->Corresponding_Pseudogene) {
        push (@notes,'This gene is thought to be a pseudogene');
    }

    if ($object->CGC_name || $object->Other_name) {
        if (my @contained_in = $object->In_cluster) {
#####      my $cluster = join ' ',map{a({-href=>Url('gene'=>"name=$_")},$_)} @contained_in;
        my $cluster = join(' ',@contained_in);
        push @notes,"This gene is contained in gene cluster $cluster.\n";
    }

#####    push @notes,map { GetEvidence(-obj=>$_,-dont_link=>1) } $object->Remark if $object->Remark;
    push @notes,$object->Remark if $object->Remark;
    }

    # Add a brief remark for Transposon CDS entries
    push @notes,
    'This gene is believed to represent the remnant of a transposon which is no longer functional'
    if (eval {$object->Corresponding_CDS->Method eq 'Transposon_CDS'});

    foreach (@notes) {
        $_ = ucfirst($_);
        $_ .= '.' unless /\.$/;
    }
    return \@notes;
}

sub parse_year { # pending deletion
    my $date = shift;
    $date =~ /.*(\d\d\d\d).*/;
    my $year = $1 || $date;
    return $year;
}


sub _pattern_thumbnail {
    my ($self,$ep) = @_;
    return '' unless $self->_is_cached($ep->name);
    my $terms = join ', ', map {$_->Term} $ep->Anatomy_term;
    $terms ||= "No adult terms in the database";
    return ([$ep,$terms]);
}

# Meh. This is a view component and doesn't belong here.
sub _is_cached {
    my ($self,$ep) = @_;
    my $WORMVIEW_IMG = '/usr/local/wormbase/html/images/expression/assembled/';
    return -e $WORMVIEW_IMG . "$ep.png";
}



sub _y2h_data { # pending deletion
    my ($self,$object,$limit,$c) = @_;
    my %tags = ('YH_bait'   => 'Target_overlapping_CDS',
                'YH_target' => 'Bait_overlapping_CDS');

    my %results;
    foreach my $tag (keys %tags) {
        if (my @data = $object->$tag) {

    # Map baits/targets to CDSs
            my $subtag = $tags{$tag};
            my %seen = ();
            foreach (@data) {
            my @cds = $_->$subtag;

            unless (@cds) {
                my $try_again = ($subtag eq 'Bait_overlapping_CDS') ? 'Target_overlapping_CDS' : 'Bait_overlapping_CDS';
                @cds = $_->$try_again;
            }

            unless (@cds) {
                my $try_again = ($subtag eq 'Bait_overlapping_CDS') ? 'Bait_overlapping_gene' : 'Target_overlapping_gene';
                my $new_gene = $_->$try_again;
                @cds = $new_gene->Corresponding_CDS if $new_gene;
            }

            foreach my $cds (@cds) {
                push @{$seen{$cds}},$_;
            }
            }

            my $count = 0;
            for my $cds (keys %seen){
                my ($y2h_ref,$count);
                my $str = "See: ";
                for my $y2h (@{$seen{$cds}}) {
                    $count++;
                    # If we are limiting for the main page, append a link to "more"
                    last if ($limit && $count > $limit);
                    #      $str    .= " " . $c->object2link($y2h);
                    $str    .= " " . $y2h;
                    $y2h_ref  = $y2h->Reference;
                }
                if ($limit && $count > $limit) {
                #      my $link = DisplayMoreLink(\@data,'y2h',undef,'more',1);
                #      $link =~ s/[\[\]]//g;
                #      $str .= " $link";
                }
                my $dbh = $self->service('acedb');
                my $k_cds = $dbh->fetch(CDS => $cds);
                #    push @{$results{$tag}}, [$c->object2link($k_cds) . " [" . $str ."]", $y2h_ref];
                push @{$results{$tag}}, [$k_cds . " [" . $str ."]", $y2h_ref];
            }
        }
    }
    return (\@{$results{'YH_bait'}},\@{$results{'YH_target'}});
}


=pod
# This is one big ugly hack job, evidence is handled by _get_evidence in API/Object.pm
sub _go_evidence_code { # pending deletion
    my ($self,$term) = @_;
    my @type      = $term->col;
    my @evidence  = $term->right->col if $term->right;
    my @results;
    foreach my $type (@type) {
    my $evidence = '';

    for my $ev (@evidence) {
        my $desc;
        my (@supporting_data) = $ev->col;

        # For IMP, this is semi-formatted text remark
        if ($type eq 'IMP' && $type->right eq 'Inferred_automatically') {
        my (%phenes,%rnai);
        foreach (@supporting_data) {
            my @row;
            $_ =~ /(.*) \(WBPhenotype(.*)\|WBRNAi(.*)\)/;
            my ($phene,$wb_phene,$wb_rnai) = ($1,$2,$3);
            $rnai{$wb_rnai}++ if $wb_rnai;
            $phenes{$wb_phene}++ if $wb_phene;
        }
#    $evidence .= 'via Phenotype: '
#      #          . join(', ',map { a({-href=>ObjectLink('phenotype',"WBPhenotype$_")},$_) }
#      . join(', ',map { a({-href=>Object2URL("WBPhenotype$_",'phenotype')},$_) }
#
#         keys %phenes) if keys %phenes > 0;

        $evidence .= 'via Phenotype: '
            . join(', ',         keys %phenes) if keys %phenes > 0;

        $evidence .= '; ' if $evidence && keys %rnai > 0;

#    $evidence .= 'via RNAi: '
#      . join(', ',map { a({-href=>Object2URL("WBRNAi$_",'rnai')},$_) }
#         keys %rnai) if keys %rnai > 0;
        $evidence .= 'via RNAi: '
            . join(', ', keys %rnai) if keys %rnai > 0;

        next;
        }

        my @seen;

        foreach (@supporting_data) {
        if ($_->class eq 'Paper') {  # a paper
#      push @seen,ObjectLink($_,build_citation(-paper=>$_,-format=>'short'));

            push @seen,$_;
        } elsif ($_->class eq 'Person') {
            #          push @seen,ObjectLink($_,$_->Standard_name);
            next;
        } elsif ($_->class eq 'Text' && $ev =~ /Protein/) {  # a protein
#      push @seen,a({-href=>sprintf(Configuration->Protein_links->{NCBI},$_),-target=>'_blank'},$_);
        } else {
#      push @seen,ObjectLink($_);
            push @seen,$_;
        }
        }
        if (@seen) {
        $evidence .= ($evidence ? ' and ' : '') . "via $desc ";
        $evidence .= join('; ',@seen);
        }
    }


    # Return an array of arrays, containing the go evidence code (IMP, IEA) and its source (RNAi, paper, curator, etc)
    push @results,[$type,($type eq 'IEA') ? 'via InterPro' : $evidence];
    }
    #my @proteins = $term->at('Protein_id_evidence');
    return @results;
}

=cut

sub _build_hash {
    open my $fh, '<', $_[0] or die $!;

    return { map { chomp; split /=>/, $_, 2 } <$fh> };
}

# helper method, retrieve public name from objects
sub _public_name {

    my ($self,$object) = @_;
    my $common_name;
    my $class = eval{$object->class} || "";

    if ($class =~ /gene/i) {
        $common_name =
        $object->Public_name
        || $object->CGC_name
        || $object->Molecular_name
        || eval { $object->Corresponding_CDS->Corresponding_protein }
        || $object;
    }
    elsif ($class =~ /protein/i) {
        $common_name =
        $object->Gene_name
        || eval { $object->Corresponding_CDS->Corresponding_protein }
        ||$object;
    }
    else {
        $common_name = $object;
    }

    my $data = $common_name;
    return "$data";


}

#######################################
#
# OBSOLETE METHODS?
#
#######################################

# Fetch all proteins associated with a gene.
## NB: figure out the naming convention for proteins

# NOTE: this method is not used
# sub proteins {
#     my $self   = shift;
#     my $object = $self->object;
#     my $desc = 'proteins related to gene';

#     my @cds    = $object->Corresponding_CDS;
#     my @proteins  = map { $_->Corresponding_protein } @cds;
#     @proteins = map {$self->_pack_obj($_, $self->public_name($_, $_->class))} @proteins;

#     return { description => 'proteins encoded by this gene',
#          data        => \@proteins };
# }


# # Fetch all CDSs associated with a gene.
# ## figure out naming convention for CDs

# # NOTE: this method is not used
# sub cds {
#     my $self   = shift;
#     my $object = $self->object;
#     my @cds    = $object->Corresponding_CDS;
#     my $data_pack = $self->basic_package(\@cds);

#     return { description => 'CDSs encoded by this gene',
#          data        => $data_pack };
# }



# # Fetch Homology Group Objects for this gene.
# # Each is associated with a protein and we should probably
# # retain that relationship

# # NOTE: this method is not used
# # TH: NOT YET CLEANED UP
# sub kogs {
#     my $self   = shift;
#     my $object = $self->object;
#     my @cds    = $object->Corresponding_CDS;
#     my %data;
#     my %data_pack;

#     if (@cds) {
#     my @proteins  = map {$_->Corresponding_protein(-fill=>1)} @cds;
#     if (@proteins) {
#         my %seen;
#         my @kogs = grep {$_->Group_type ne 'InParanoid_group' } grep {!$seen{$_}++}
#              map {$_->Homology_group} @proteins;
#         if (@kogs) {

#             $data_pack{$object} = \@kogs;
#             $data{'data'} = \%data_pack;

#         } else {

#             $data_pack{$object} = 1;

#         }
#     }
#     } else {
#         $data_pack{$object} = 1;
#     }

#     $data{'description'} = "KOGs related to gene";
#      return \%data;
# }

__PACKAGE__->meta->make_immutable;

1;
