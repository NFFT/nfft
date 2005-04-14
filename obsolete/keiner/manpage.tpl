{(   -*- nroff -*-

autogen5 template script=updater

)}
{(#

1999-06-26 original version, Jim Van Zandt <jrv@vanzandt.mv.com>
1999-06-28 modified to be template  Bruce Korb <autogen@linuxbox.com>
1999-07-04 nroff tweaks, Jim Van Zandt <jrv@vanzandt.mv.com>
2001-08-25 converted to autogen5 template, Jim Van Zandt <jrv@vanzandt.mv.com>

)}/^.\\" @synopsis@$/,/^.\\" @@$/{
  if (!synopsis_updated) {
    print ".\\\" @synopsis@\n\
.B {(
  IF (exist? "prog_file_name") )}{(prog_file_name)}{(
  ELSE)}{(prog_name)}{(
  ENDIF)}{(

  IF (exist? "flag.value") )}{(
    IF (exist? "long_opts") )}{(

      # * * * * * * * * * * * * * * * * * * * * * * * * *
      #
      #  Mixture of short (flag) options and long options
      #
      )}\n\
.\\\" Mixture of short (flag) options and long options\n\
.RB [ -\\fIflag\\fP \" [\\fIvalue\\fP]]... [\" --\\fIopt-name\\fP [ = \"| ]\\fIvalue\\fP]...\"{(

    ELSE )}{(

      # * * * * * * * * * * * * * * * * * * * * * * * * *
      #
      #  Short (flag) options only
      #
      )}
.\\\" Short (flag) options only\n\
.RB [ -\\fIflag\\fP \" [\\fIvalue\\fP]]...\"{(
    ENDIF )}{(

  ELIF (exist? "long_opts") )}{(

      # * * * * * * * * * * * * * * * * * * * * * * * * *
      #
      #  Long options only
      #
      )}\n\
.\\\" Long options only\n\
.RB [ --\\fIopt-name\\fP [ = \"| ] \\fIvalue\\fP]...\"{(

  ELIF (not (exist? "argument") )}{(

      # * * * * * * * * * * * * * * * * * * * * * * * * *
      #
      #  All arguments are named options.
      #
      )}
.\\\" All arguments are named options.\n\
.RI [ opt-name \"[\\fB=\\fP\" value ]]...\n\
.PP\n\
All arguments are named options.{(

  ELSE )}{(
    (error "Named option programs cannot have arguments") )}{(
  ENDIF )}{(

  IF (exist? "argument") )}
.br\n\
.in +8\n\
{(argument)}{(
  ENDIF )}";
  print ".\\\" @@"
  }
  synopsis_updated=1;
  next;
}
/^.\\" @options@$/,/^.\\" @@$/{
  if (!options_updated) {
    print ".\\\" @options@";
{(

#  When there are no flag characters and long options are disallowed,
   the program options are all named options and they are handled
   without the marker character (hyphen). )}{(

FOR flag

)}{(
  #  Skip the documentation options!
  #
  )}{(
  IF (not (exist? "documentation")) )}print "\
.TP{(
      IF (exist? "value") )}{(
        IF (exist? "long_opts") )}{(

          # * * * * * * * * * * * * * * * * * * * *
          *
          *  The option has a flag value (character) AND
          *  the program uses long options
          *
          )}\n\
.BR -{(value)}{(
          IF (not (exist? "arg_name")) )} \", \" --{(
          ELSE )} \" \\fI{(arg_name)}\\fP, \" --{(
          ENDIF )}{(% name (string-tr! "%s" "A-Z_^" "a-z--") )}{(
          IF (exist? "arg_name") )}{(
            IF ( = (len "flag_arg") 0 )} [ ={(
            ELSE )} \" \" {(
            ENDIF )}\\fI{(arg_name)}\\fP{(
            IF ( = (len "flag_arg") 0 )} ]{(
            ENDIF )}{(
          ENDIF)}{(


        ELSE   )}{(

          # * * * * * * * * * * * * * * * * * * * *
          *
          *  The option has a flag value (character) BUT
          *  the program does _NOT_ use long options
          *
          )}\n\
.BR -{(value)}{(
          IF (exist" "arg_name) )}{(
            IF ( = (len "flag_arg") 0 )}[{(
            ENDIF )} \"\\fI{(arg_name)}\\fP{(
            IF ( = (len "flag_arg") 0 )}\"]\"{(
            ENDIF )}{(
          ENDIF)}{(
        ENDIF  )}{(


      ELSE  value does not exist -- named option only  )}{(

        IF (not (exist? "long_opts") )}{(

          # * * * * * * * * * * * * * * * * * * * *
          *
          *  The option does not have a flag value (character).
          *  The program does _NOT_ use long options either.
          *  Special magic:  All arguments are named options.
          *
          )}\n\
.BR {((string-downcase (string-tr! "name" "#_^" "#--")))}{(
        IF arg_name _exist)} {(
          IF ( = (len "flag_arg") 0) )} [ ={(
          ELSE )} \"\" {(
          ENDIF )}\\fI{(arg_name)}\\fP{(
          IF (= (len "flag_arg") 0) )}]{(
          ENDIF )}{(
        ENDIF)}{(


        ELSE   )}{(

          # * * * * * * * * * * * * * * * * * * * *
          *
          *  The option does not have a flag value (character).
          *  The program, instead, only accepts long options.
          *
          )}\n\
.BR --{((string-downcase (string-tr! "name" "#_^" "#--")))}{(
          IF (exist? "arg_name"))} \"{(
            IF ( = (len "flag_arg") 0 ) )}[{(
            ENDIF )}\" \"=\\fI{(arg_name)}\\fP{(
            IF ( = (len "flag_arg") 0 ) )}]{(
            ENDIF )}\"{(
          ENDIF)}{(
        ENDIF  )}{(
      ENDIF )}\n\
{(descrip)}{(
      IF (or (exist? "min")
	     (exist? "max")
	     (exist? "enabled")
	     (exist? "no_preset")
	     (exist? "setable")
	     (exist? "equivalence")
	     (exist? "flags_must")
	     (exist? "flags_cant")) )}.\n\
Special attributes apply to this option.\n\
See the invocation section in the info document.{(
      ENDIF)}\n\
{(    IF (exist? "doc"))}.sp\n\
{(    ENDIF)}{( 
# SED
   -e convert @code into \fB...\fP phrases
   -e convert @file into \fI...\fP phrases
   -e convert @var into \fB...\fP phrases
   -e Remove the '@' prefix from curly braces
   -e Indent example regions
   -e Delete example command
   -e Replace "end example" command with ".br"
   -e Replace "@*" command with ".br"
   -e Add "\n\" to the end of each line.

   NB:  backslashes are interpreted four times: by AutoGen, sed, awk 
	and nroff. Thus, for nroff to see one backslash, use eight (8)!!

   ))}{(
doc (shell (printf "sed					\
 -e's;@code{\\([^}]*\\)\};\\\\\\\\fB\\1\\\\\\\\fP;g'	\
 -e's;@file{\\([^}]*\\)\};\\\\\\\\fI\\1\\\\\\\\fP;g'	\
 -e's;@var{\\([^}]*\\)\};\\\\\\\\fB\\1\\\\\\\\fP;g'	\
 -e's/@\\([{}]\\)/\\1/'					\
 -e'/@ *example/,/@ *end *example/s/^/    /'		\
 -e'/^ *@ *example/d'					\
 -e's/^ *@ *end *example/.br/'				\
 -e's/^@\\*/.br/'					\
 -e's/$/\\\\n\\\\/' <<'_EOF_'\n%s\n_EOF_")) )}{(

  ENDIF (not (exist? "documentation")) )}.br";
{(

ENDFOR flag


)}print "\
.TP\n\
.BR {(IF (exist? "flag.value") )}\\-? , \" --{(
      ELIF (exist? "long_opt") )}--{(ENDIF)}help\n\
Display usage information and exit.\n\
.TP\n\
.BR {(IF (exist? "flag.value") )}-! , \" --{(
      ELIF (exist? "long_opt") )}--{(ENDIF)}more-help\n\
Extended usage information passed thru pager.{(


IF (exist? "homerc") )}\n\
.TP\n\
.BR {(IF (exist? "flag.value") )}-> \" \\fIrcfile\\fP, --\" {(
      ELIF (exist? "long_opt") )}--{(ENDIF)}save-opts \"[=\\fIrcfile\\fP] \n\
Save the option state to \\fIrcfile\\fP.\n\
.TP\n\
.BR {(IF (exist? "flag.value") )}-< \" \\fIrcfile\\fP, --\" {(
      ELIF (exist? "long_opt") )}--{(ENDIF)}load-opts \"=\\fIrcfile\\fP\"\n\
Load options from \\fIrcfile\\fP.\n\
.TP\n\
.BR  --no-load-opts\n\
Disable loading options from an rc file.{(
ENDIF)}{(


IF (exist? "version") )}\n\
.TP\n\
.BR {(IF (exist? "flag.value") )}\\-v \" [\" v | c | n \"], \" --{(
      ELIF (exist? "long_opt") )}--{(ENDIF)}version [ =v | c | n ]\n\
Output version of program and exit.  The default mode is `v', a simple\n\
version.  The `c' mode will print copyright information and `n' will\n\
print the full copyright notice.{(
ENDIF


)}";
  print ".\\\" @@"
  }
  options_updated=1;
  next;
}
{print}
