package WormBase::Datomic;

use strict;

use edn;
use HTTP::Tiny;

# Some boiler-plate for querying the Datomic REST service.
# Currently also holds default URI and DB-alias until we
# Think of a better place to keep them!

sub new {
    my ($class, $uri, $db) = @_;
    $uri = $uri || 'http://localhost:4664/api/query';
    $db = $db || 'ace/wb244-imp2';
    bless {
	uri   => $uri,
	db    => $db,
        http  => HTTP::Tiny->new
    }, $class;
}

sub query {
    my $self = shift;
    my $query = shift;
    my $params = edn::write([
	{'db/alias' => $self->{'db'}},
	@_
    ]);
    my $resp = $self->{'http'}->post_form(
	$self->{'uri'},
	{q     => $query,
         args  => $params}
    );
    die $resp->{'content'} unless $resp->{'status'} == 200;
    edn::read($resp->{'content'});
}

1;
