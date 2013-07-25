
package WormBase::API::Role::Interaction;

use Moose::Role;
use Data::Dumper;

#######################################################
#
# Attributes
#
#######################################################

has '_interactions' => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build__interactions',
);

sub _build__interactions {
    my $self = shift;
	my $data;
    return $self->_get_interactions($data, 0);
}

#######################################
#
# The Interactions Widget
#   template: classes/gene/interactions.tt2
#
#######################################

=head2 Interactions

=cut

=head3 interactions

This method returns a data structure containing the 
a data table of gene and protein interactions. Ask us
to increase the granularity of this method!

=over

=item PERL API

 $data = $model->interactions();

=item REST API

B<Request Method>

GET

B<Requires Authentication>

No

B<Parameters>

An object ID (WBGene00006763)

B<Returns>

=over 4

=item *

200 OK and JSON, HTML, or XML

=item *

404 Not Found

=back

B<Request example>

curl -H content-type:application/json http://api.wormbase.org/rest/field/[class]/[object]/interactions

B<Response example>

=cut

sub interactions  {
    my $self   = shift;
    my $object = $self->object;
    my $class = $object->class;

    my @edges = values %{$self->_interactions->{edgeVals}};

    my $results = $self->_get_interactions($self->_interactions, 1, 1);
    my @edges_all = values %{$results->{edgeVals}};

    return {
        description => 'genetic and predicted interactions',
        data        => $results->{showall} ? {    
                            edges => @edges ? \@edges : undef,
                            types => $results->{types},
                            nodes => $results->{nodes},
                            showall => $results->{showall},
                            ntypes => $results->{ntypes},
                            edges_all => @edges_all ? \@edges_all : undef,
                            class => $class
                       } : { edges => \@edges },
    };

}

=head3 interaction_details

This method takes the interaction network built on first load and
searches for additional interactions for the same set of nodes.
Returns a data structure in the same format as the one taken in.

=over

=item PERL API

 $data = $model->interaction_details();

=item REST API

B<Request Method>

GET

B<Requires Authentication>

No

B<Parameters>

An object ID (WBGene00006763)

B<Returns>

=over 4

=item *

200 OK and JSON, HTML, or XML

=item *

404 Not Found

=back

B<Request example>

curl -H content-type:application/json http://api.wormbase.org/rest/field/[class]/[object]/interaction_details

B<Response example>

=cut

sub interaction_details {
    my $self = shift;
    my $results = $self->_get_interactions($self->_interactions, 1);

    my @edges = values %{$results->{edgeVals}};

    return {
		description	=> 'addtional nearby interactions',
		data		=> {    edges => @edges ? \@edges : undef,
                            types => $results->{types},
                            nodes => $results->{nodes},
                            ntypes => $results->{ntypes},
                            showall => $results->{showall},
                       },
    };
}

############################################################
#
# Private Methods
#
############################################################

sub _get_interactions {
    my ($self, $data, $nearby, $from_table) = @_;
    my $object = $self->object;
    my @objects;
    #determine object type and extract interactions accordingly
    if ($nearby){ 
      @objects = map {$_->Interaction} grep {($_->class =~ /gene/i) && ($data->{nodes}{"$_"}{predicted} == 0)} values %{$data->{nodes_obj}} 
    } elsif ($object->class =~ /gene|wbprocess/i ) { # ELSEIF PROCESS, COPY GENE FUNC
      @objects = $object->Interaction 
    } elsif ($object->class =~ /interaction/i ) { 
      @objects = ($object) 
    } 
    if($nearby && (scalar @objects > 3000)){
        $data->{showall} = 0;
        return $data;
    }

    $self->log->debug("nearby: $nearby, size: ", scalar @objects);

    foreach my $interaction ( @objects ) {
      next if($data->{ids}{"$interaction"});
      if ($nearby) { next if scalar grep {!defined $data->{nodes_obj}->{$_} || ($data->{nodes}{"$_"}{predicted} == 1)} map {$_->col} $interaction->Interactor; }
      my $edgeList = $self->_get_interaction_info($interaction, $nearby);
      foreach my $key (keys %{$edgeList}) {
          my ($type, $effector, $affected, $direction, $phenotype)= @{$edgeList->{$key}};
          next unless($effector);
          next unless ($effector->class =~ /Molecule|Gene|Rearrangement|Text/ && $affected->class =~ /Molecule|Gene|Rearrangement|Text/);

          $data->{nodes}{"$effector"} ||= $self->_pack_obj($effector);
          $data->{nodes}{"$affected"} ||= $self->_pack_obj($affected);

          $data->{nodes}{"$effector"}{predicted} = ($type eq 'Predicted') ? $data->{nodes}{"$effector"}{predicted} // 1 : 0;
          $data->{nodes}{"$affected"}{predicted} = ($type eq 'Predicted') ? $data->{nodes}{"$affected"}{predicted} // 1 : 0;

          $data->{nodes_obj}{"$effector"} = $effector;
          $data->{nodes_obj}{"$affected"} = $affected;
          $data->{ids}{"$interaction"}=1;
          $data->{types}{"$type"}=1;
            
          my $ntype1 = $data->{nodes}{"$effector"}->{class};
          my $ntype2 = $data->{nodes}{"$affected"}->{class};
          $data->{ntypes}{"$ntype1"}=1;
          $data->{ntypes}{"$ntype2"}=1;

          my $phenObj = $self->_pack_obj($phenotype);
          my $key = "$effector $affected $type";
          my $key2 = "$affected $effector $type";

          if ($phenotype) {
            $data->{phenotypes}{"$phenotype"} ||= $phenObj;
            $key .= " $phenotype" if $phenotype;
            $key2 .= " $phenotype" if $phenotype;
          }
          
          my $packInteraction = $self->_pack_obj($interaction);
          my @papers = map { $self->_pack_obj($_) } $interaction->Paper;
          
          if (exists $data->{edgeVals}{$key}){
            push @{$data->{edgeVals}{$key}{interactions}}, $packInteraction;
            push @{$data->{edgeVals}{$key}{citations}}, @papers;
          } elsif (exists $data->{edgeVals}{$key2}){
            push @{$data->{edgeVals}{$key2}{interactions}}, $packInteraction;
            push @{$data->{edgeVals}{$key2}{citations}}, @papers;
          } else {
            my @interacArr = ($packInteraction);
            $data->{edgeVals}{$key} = {
                interactions=> @interacArr ? \@interacArr : undef,
                citations	=> @papers ? \@papers : undef,
                type	=> "$type",
                effector	=> $data->{nodes}{"$effector"},
                affected	=> $data->{nodes}{"$affected"},
                direction	=> $direction,
                phenotype	=> $phenObj,
                nearby	=> $nearby,
            };
          }
      }
    }
    $data->{showall} = scalar keys %{$data->{edgeVals}} < 100 || $nearby;
    return $data;
}

sub _get_interaction_info {
    my ($self, $interaction, $nearby) = @_;
    my %results;
    my $object = $self->object;
    my $type = $interaction->Interaction_type;
    $type = $type->right ? $type->right . '' : "$type";
    return \%results if(($type eq 'No_interaction'));
    $type =~ s/_/ /g;
    if ($type eq 'Regulatory') {
    if ( my $reg_result = $interaction->Regulation_result ) {
        if ("$reg_result" =~ /^(.*tive)_regulate$/) { $type = $1 . "ly Regulates" }
        elsif ("$reg_result" eq 'Does_not_regulate') { $type = "Does Not Regulate" }
    }
    }
    # Filter low confidence predicted interactions.
    return \%results if ($interaction->Log_likelihood_score || 1000) <= 1.5 && $object->class ne 'Interaction' && $type =~ m/predicted/i;

    my $phenotype = $interaction->Interaction_phenotype;
    my ( @effectors, @affected, @others );
    my %elements;
    foreach my $intertype ($interaction->Interactor) {
      my $count = 0;
      foreach my $interactor ($intertype->col) {
          my @interactors = $intertype->col;
          my @tags = eval { $interactors[$count++]->col };

          my %info;
          $info{obj} = $interactor;
          if ( @tags ) {
            map { $info{"$_"} = $_->at; } @tags;
            if ("$intertype" eq 'Interactor_overlapping_gene') {
                my $role = $info{Interactor_type};
                if ($role && $role =~ /Effector|.*regulator/) {   push @effectors, $interactor }
                elsif ($role && $role =~ /Affected|.*regulated/)  { push @affected, $interactor }
                else { push @others, $interactor }
            } else {
                my $corresponding_gene = $self->_get_gene($interactor, "$intertype");
                if ($corresponding_gene) { @{$results{"$interactor $corresponding_gene"}} = ('Associated Product', $interactor, $corresponding_gene, 'Other') } 
                else { push @others, $interactor }
            }
          } else { push @others, $interactor }
      }
    }
    if (@effectors || @affected) {
      foreach my $obj (@effectors, @others) {
          foreach my $obj2 (@affected) {
            next if "$obj" eq "$obj2";
            if (!$nearby && $object->class ne 'Interaction' && $object->class ne 'WBProcess') { 
              next unless ("$obj" eq "$object" || "$obj2" eq "$object")
            };
            @{$results{"$obj $obj2"}} = ($type, $obj, $obj2, 'Effector->Affected', $phenotype);
          }
      }
    } else {
      foreach my $obj (@others) {
          foreach my $obj2 (@others) {
            next if "$obj" eq "$obj2";
            if (!$nearby && $object->class ne 'Interaction') { 
              next unless ("$obj" eq "$object" || "$obj2" eq "$object")
            };
            my @objs = ("$obj", "$obj2");
            my $str = join(' ', sort @objs); 
            @{$results{"$str"}} = ($type, $obj, $obj2, 'non-directional', $phenotype);
          }
      }
    }
    return \%results;
}

sub _get_gene {
    my ($self, $obj, $type) = @_;
    if ($type eq 'Interactor_overlapping_CDS') { return $obj->Gene }
    elsif ($type eq 'Interactor_overlapping_protein') { return $obj->Corresponding_CDS->Gene if $obj->Corresponding_CDS }
    elsif ($type eq 'PCR_interactor') {
      my $corr_gene = $obj->Overlaps_CDS->Gene if $obj->Overlaps_CDS;
      $corr_gene ||= $obj->Overlaps_transcript->Gene if $obj->Overlaps_transcript;
      $corr_gene ||= $obj->Overlaps_pseudogene->Gene if $obj->Overlaps_pseudogene;
      return $corr_gene;
    } elsif ($type eq 'Sequence_interactor') {
      my $corr_gene = $obj->Matching_CDS->Gene if $obj->Matching_CDS;
      $corr_gene ||= $obj->Matching_transcript->Gene if $obj->Matching_transcript;
      $corr_gene ||= $obj->Matching_pseudogene->Gene if $obj->Matching_pseudogene;
      return $corr_gene;
    }
}

1;
