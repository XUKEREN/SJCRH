package XMLTools;

use strict;
use Configurable;
use XML::Parser;

@XMLTools::ISA = qw(Configurable Exporter);
@XMLTools::EXPORT_OK = qw();

use MethodMaker qw(
		   file
		   tree
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);

  if (my $fn = $self->file) {
    die "can't find file $fn" unless -s $fn;
    my $xp = new XML::Parser(Style => "Tree");
    $self->tree($xp->parsefile($fn));
  }

  return $self;
}

sub find_tags {
  #
  # traverse an XML tree, returning nodes for a given tag
  #
  my ($self, %options) = @_;
  my $search_tag = $options{"-tag"} || die "specify -tag";
  $search_tag = lc($search_tag);
  my $tree = $options{"-tree"} || $self->tree();
  my $results = $options{"-results"} || [];
  
  my ($tag, $content) = @{$tree};
  if (lc($tag) eq $search_tag) {
    push @{$results}, $tree;
  } elsif ($tag eq "0") {
    # pseudotag for plain text (ignore)
  } else {
    my ($attrs, @pairs) = @{$content};
    my $len = length @pairs;
    for (my $i=0; $i < @pairs; $i += 2) {
      $self->find_tags(
		       "-tag" => $search_tag,
		       "-tree" => [ @pairs[$i, $i+1] ],
		       "-results" => $results
		      );
      # recurse
    }
  }
  return $results;
}

sub find_tag_text {
  #
  #  traverse a tree, returning text content for a tag
  #
  my ($self, %options) = @_;
  my $search_tag = $options{"-tag"} || die "specify -tag";
  $search_tag = lc($search_tag);
  my $tree = $options{"-tree"} || $self->tree();
  my $results = $options{"-results"} || [];
  my $capture = $options{"-capture"};
  my $text = $options{"-text"};
  
  my ($tag, $content) = @{$tree};
  $capture = 1 if (lc($tag) eq $search_tag);

  if ($tag eq "0") {
    # pseudotag for plain text
    push @{$results}, $content if $capture;
  } else {
    my ($attrs, @pairs) = @{$content};
    my $len = length @pairs;
    for (my $i=0; $i < @pairs; $i += 2) {
      $self->find_tag_text(
			   "-tag" => $search_tag,
			   "-tree" => [ @pairs[$i, $i+1] ],
			   "-results" => $results,
			   "-capture" => $capture
			  );
      # recurse
    }
  }

  if ($options{"-text"}) {
    return join " ", @{$results};
    # TO DO: maybe no extra whitespace if prev entry already ends with it?
  } else {
    return $results;
  }
}

sub get_node_attr_hash {
  my ($self, $node) = @_;
  return $node->[1]->[0];
}

sub get_node_attribute {
  my ($self, $node, $attribute) = @_;
  return $node->[1]->[0]->{$attribute};
}

sub get_tag_attribute {
  my ($self, %options) = @_;
  my $nodes = $self->find_tags(%options);
  die "nodes found != 1" unless @{$nodes} == 1;
  return $self->get_node_attribute($nodes->[0], $options{"-attr"});
}


1;
