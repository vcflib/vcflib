use strict;
use warnings;

package Local::vcflib::Test;
use base 'Exporter';

use File::Basename qw< dirname >;
use IPC::Run3 qw< run3 >;
use Test::More;

our @EXPORT = qw( run run_ok );
our $BIN    = dirname(__FILE__) . "/../../../../bin";

sub run {
    my ($run, $stdin)    = @_;
    my ($command, @opts) = @$run;
    run3(["$BIN/$command", @opts], \$stdin, \(my $stdout), \(my $stderr));
    return ($stdout, $stderr, $?);
}

sub run_ok {
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my ($stdout, $stderr, $exit) = run(@_);
    ok $exit >> 8 == 0, "exit code"
        or diag "error running command: " . join(" ", @{$_[0]}) . "\n"
               ."with input:\n$_[1]\n--\n"
               ."exit code = " . ($exit >> 8) . " (system() return value = $exit)\n"
               ."stderr = \n$stderr";
    return ($stdout, $stderr);
}

1;
