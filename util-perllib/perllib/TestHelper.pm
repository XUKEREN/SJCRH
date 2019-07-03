package TestHelper; 

#use base 'Test::Builder::Module'; 
use Test::More; 

my $CLASS = __PACKAGE__; 

sub new_ok_hash { 
    my $tb = Test::More->builder;
    $tb->croak("new_ok_hash() must be given at least a class") unless @_;

    my( $class, $args, $object_name ) = @_;
 
    $args ||= {};
    $object_name = "The object" unless defined $object_name;

    my $obj;
    my( $success, $error ) = $tb->_try( sub { $obj = $class->new(%$args); 1 } );
    if($success) {
        local $Test::Builder::Level = $Test::Builder::Level + 1;
        &Test::More::isa_ok($obj, $class, $object_name);
    }
    else {
        $tb->ok( 0, "new() died" );
        $tb->diag("    Error was:  $error");
    }

    return $obj;

}

sub ok_or_die{
	my ($test, $name) = @_; 
	my $tb = Test::More->builder; 
	
	#print STDERR "$test, $name\n"; 
	my $res = $tb->ok( $test, $name ); 
	if ($res != 1){
		#rather than die() let's use the testing framework's ability to bail 
		#die "Critical test failed";
		BAIL_OUT("Critical test failed");  
	}
	
	return $res; 	

}
