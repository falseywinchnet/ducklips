// Pronunciation: a hand-written lexicon of frequent and irregular words,
// then letter-to-sound rules in the style of Elovitz, Johnson, McHugh and
// Shore (NRL report 7948, 1976), then a stress heuristic.

#include "voice/text.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>

namespace voice {

namespace {

// word PHONES (stress digits on vowels)
constexpr const char* kLexicon = R"(a AH0
about AH0 B AW1 T
above AH0 B AH1 V
after AE1 F T ER0
again AH0 G EH1 N
against AH0 G EH1 N S T
ago AH0 G OW1
air EH1 R
alfa AE1 L F AH0
all AO1 L
allow AH0 L AW1
almost AO1 L M OW2 S T
alone AH0 L OW1 N
along AH0 L AO1 NG
alpha AE1 L F AH0
already AO0 L R EH1 D IY0
also AO1 L S OW0
always AO1 L W EY0 Z
am AE1 M
among AH0 M AH1 NG
an AE1 N
and AE1 N D
another AH0 N AH1 DH ER0
any EH1 N IY0
anything EH1 N IY0 TH IH2 NG
are AA1 R
around ER0 AW1 N D
as AE1 Z
ask AE1 S K
at AE1 T
away AH0 W EY1
back B AE1 K
be B IY1
because B IH0 K AH1 Z
been B IH1 N
before B IH0 F AO1 R
began B IH0 G AE1 N
being B IY1 IH0 NG
believe B IH0 L IY1 V
below B IH0 L OW1
between B IH0 T W IY1 N
both B OW1 TH
boy B OY1
bravo B R AA1 V OW0
bread B R EH1 D
break B R EY1 K
breath B R EH1 TH
brother B R AH1 DH ER0
brown B R AW1 N
build B IH1 L D
busy B IH1 Z IY0
but B AH1 T
buy B AY1
by B AY1
call K AO1 L
came K EY1 M
can K AE1 N
charlie CH AA1 R L IY0
child CH AY1 L D
children CH IH1 L D R AH0 N
city S IH1 T IY0
clown K L AW1 N
come K AH1 M
could K UH1 D
country K AH1 N T R IY0
cow K AW1
crowd K R AW1 D
crown K R AW1 N
day D EY1
delta D EH1 L T AH0
did D IH1 D
do D UW1
does D AH1 Z
done D AH1 N
door D AO1 R
down D AW1 N
drown D R AW1 N
each IY1 CH
early ER1 L IY0
earth ER1 TH
echo EH1 K OW0
eight EY1 T
eighteen EY0 T IY1 N
eighty EY1 T IY0
eleven IH0 L EH1 V AH0 N
enough IH0 N AH1 F
even IY1 V IH0 N
evening IY1 V N IH0 NG
ever EH1 V ER0
every EH1 V R IY0
everything EH1 V R IY0 TH IH2 NG
eye AY1
eyes AY1 Z
father F AA1 DH ER0
few F Y UW1
fifteen F IH0 F T IY1 N
fifty F IH1 F T IY0
find F AY1 N D
first F ER1 S T
five F AY1 V
flower F L AW1 ER0
for F AO1 R
forty F AO1 R T IY0
four F AO1 R
fourteen F AO0 R T IY1 N
foxtrot F AA1 K S T R AA2 T
friend F R EH1 N D
from F R AH1 M
frown F R AW1 N
full F UH1 L
give G IH1 V
go G OW1
goes G OW1 Z
golf G AO1 L F
gone G AO1 N
good G UH1 D
gown G AW1 N
great G R EY1 T
had HH AE1 D
half HH AE1 F
has HH AE1 Z
have HH AE1 V
he HH IY1
head HH EH1 D
heard HH ER1 D
heart HH AA1 R T
hello HH AH0 L OW1
her HH ER1
here HH IY1 R
him HH IH1 M
his HH IH1 Z
home HH OW1 M
hotel HH OW0 T EH1 L
hour AW1 ER0
house HH AW1 S
how HH AW1
howl HH AW1 L
human HH Y UW1 M AH0 N
hundred HH AH1 N D R AH0 D
i AY1
if IH1 F
in IH1 N
india IH1 N D IY0 AH0
into IH1 N T UW0
is IH1 Z
it IH1 T
its IH1 T S
juliet JH UW1 L IY0 EH2 T
juliett JH UW1 L IY0 EH2 T
just JH AH1 S T
kilo K IY1 L OW0
knew N UW1
know N OW1
known N OW1 N
language L AE1 NG G W IH0 JH
large L AA1 R JH
last L AE1 S T
learn L ER1 N
leave L IY1 V
light L AY1 T
like L AY1 K
lima L IY1 M AH0
listen L IH1 S AH0 N
little L IH1 T AH0 L
live L IH1 V
long L AO1 NG
look L UH1 K
love L AH1 V
made M EY1 D
make M EY1 K
mama M AA1 M AH0
man M AE1 N
many M EH1 N IY0
may M EY1
me M IY1
mean M IY1 N
men M EH1 N
might M AY1 T
mike M AY1 K
mind M AY1 N D
money M AH1 N IY0
more M AO1 R
morning M AO1 R N IH0 NG
most M OW1 S T
mother M AH1 DH ER0
mouth M AW1 TH
move M UW1 V
much M AH1 CH
music M Y UW1 Z IH0 K
must M AH1 S T
my M AY1
myself M AY0 S EH1 L F
name N EY1 M
never N EH1 V ER0
new N UW1
next N EH1 K S T
night N AY1 T
nine N AY1 N
nineteen N AY0 N T IY1 N
ninety N AY1 N T IY0
no N OW1
none N AH1 N
not N AA1 T
nothing N AH1 TH IH0 NG
november N OW0 V EH1 M B ER0
now N AW1
of AH1 V
off AO1 F
often AO1 F AH0 N
oh OW1
old OW1 L D
on AA1 N
once W AH1 N S
one W AH1 N
only OW1 N L IY0
open OW1 P AH0 N
or AO1 R
oscar AA1 S K ER0
other AH1 DH ER0
our AW1 ER0
out AW1 T
over OW1 V ER0
owl AW1 L
own OW1 N
papa P AH2 P AA1
people P IY1 P AH0 L
person P ER1 S AH0 N
place P L EY1 S
please P L IY1 Z
point P OY1 N T
power P AW1 ER0
put P UH1 T
quebec K EH2 B EH1 K
question K W EH1 S CH AH0 N
quiet K W AY1 AH0 T
read R IY1 D
really R IH1 L IY0
right R AY1 T
romeo R OW1 M IY0 OW2
said S EH1 D
same S EY1 M
saw S AO1
say S EY1
says S EH1 Z
school S K UW1 L
sea S IY1
second S EH1 K AH0 N D
see S IY1
seven S EH1 V AH0 N
seventeen S EH2 V AH0 N T IY1 N
seventy S EH1 V AH0 N T IY0
shall SH AE1 L
she SH IY1
shore SH AO1 R
should SH UH1 D
show SH OW1
shower SH AW1 ER0
sierra S IY0 EH1 R AH0
six S IH1 K S
sixteen S IH0 K S T IY1 N
sixty S IH1 K S T IY0
small S M AO1 L
so S OW1
some S AH1 M
something S AH1 M TH IH0 NG
sometimes S AH1 M T AY2 M Z
soon S UW1 N
sound S AW1 N D
speak S P IY1 K
speech S P IY1 CH
still S T IH1 L
story S T AO1 R IY0
sure SH UH1 R
take T EY1 K
talk T AO1 K
tango T AE1 NG G OW0
tell T EH1 L
ten T EH1 N
than DH AE1 N
thank TH AE1 NG K
thanks TH AE1 NG K S
that DH AE1 T
the DH AH0
their DH EH1 R
them DH EH1 M
then DH EH1 N
there DH EH1 R
these DH IY1 Z
they DH EY1
thing TH IH1 NG
things TH IH1 NG Z
think TH IH1 NG K
thirteen TH ER0 T IY1 N
thirty TH ER1 T IY0
this DH IH1 S
those DH OW1 Z
though DH OW1
thought TH AO1 T
thousand TH AW1 Z AH0 N D
three TH R IY1
through TH R UW1
time T AY1 M
to T UW1
today T AH0 D EY1
together T AH0 G EH1 DH ER0
told T OW1 L D
too T UW1
toward T AH0 W AO1 R D
towel T AW1 AX0 L
tower T AW1 ER0
town T AW1 N
twelve T W EH1 L V
twenty T W EH1 N T IY0
two T UW1
under AH1 N D ER0
uniform Y UW1 N AH0 F AO2 R M
until AH0 N T IH1 L
up AH1 P
upon AH0 P AA1 N
us AH1 S
use Y UW1 Z
very V EH1 R IY0
victor V IH1 K T ER0
voice V OY1 S
vowel V AW1 AX0 L
want W AA1 N T
was W AA1 Z
water W AO1 T ER0
way W EY1
we W IY1
well W EH1 L
went W EH1 N T
were W ER1
what W AH1 T
when W EH1 N
where W EH1 R
which W IH1 CH
while W AY1 L
whiskey W IH1 S K IY0
who HH UW1
whole HH OW1 L
whom HH UW1 M
whose HH UW1 Z
why W AY1
will W IH1 L
with W IH1 DH
without W IH0 DH AW1 T
woman W UH1 M AH0 N
women W IH1 M AH0 N
word W ER1 D
words W ER1 D Z
work W ER1 K
world W ER1 L D
would W UH1 D
write R AY1 T
x EH1 K S
xray EH1 K S R EY2
yankee Y AE1 NG K IY0
year Y IH1 R
years Y IH1 R Z
yes Y EH1 S
yet Y EH1 T
you Y UW1
young Y AH1 NG
your Y AO1 R
zero Z IH1 R OW0
zulu Z UW1 L UW0
)";

constexpr const char* kFunctionWords[] = {
    "a", "an", "the", "and", "or", "but", "nor", "of", "to", "in", "on", "at", "by", "for", "with", "from",
    "as", "into", "onto", "upon", "than", "that", "this", "these", "those", "is", "am", "are", "was", "were",
    "be", "been", "being", "do", "does", "did", "have", "has", "had", "will", "would", "shall", "should",
    "can", "could", "may", "might", "must", "it", "its", "he", "him", "his", "she", "her", "they", "them",
    "their", "we", "us", "our", "you", "your", "i", "me", "my", "if", "so", "not", "no", "there", "which",
    "who", "whom", "whose", "what", "when", "where", "while", "because", "about", "some", "any", "all",
    "just", "up", "out", "over", "then", "too", "very"};

struct Rule {
  const char* left;
  const char* match;
  const char* right;
  const char* output;
};

// Context symbols:  # one or more vowels   * one or more consonants
//   . a voiced consonant   $ one consonant followed by E or I
//   % a suffix (ER E ES ED ING ELY)   & a sibilant   @ T S R D L Z N J TH CH SH
//   ^ one consonant   + a front vowel (E I Y)   : zero or more consonants
//   space: word edge
constexpr Rule kRules[] = {
    // A
    {"", "A", " ", "AX"}, {" ", "ARE", " ", "AA R"}, {" ", "AR", "O", "AX R"}, {"", "AR", "#", "EH R"},
    {" ^", "AS", "#", "EY S"}, {"", "A", "WA", "AX"}, {"", "AW", "", "AO"}, {" :", "ANY", "", "EH N IY"},
    {"", "A", "^+#", "EY"}, {"#:", "ALLY", "", "AX L IY"}, {" ", "AL", "#", "AX L"}, {"", "AGAIN", "", "AX G EH N"},
    {"#:", "AG", "E", "IH JH"}, {"", "A", "^+:#", "AE"}, {" :", "A", "^+ ", "EY"}, {"", "A", "^%", "EY"},
    {" ", "ARR", "", "AX R"}, {"", "ARR", "", "AE R"}, {" :", "AR", " ", "AA R"}, {"", "AR", " ", "ER"},
    {"", "AR", "", "AA R"}, {"", "AIR", "", "EH R"}, {"", "AI", "", "EY"}, {"", "AY", "", "EY"},
    {"", "AU", "", "AO"}, {"#:", "AL", " ", "AX L"}, {"#:", "ALS", " ", "AX L Z"}, {"", "ALK", "", "AO K"},
    {"", "AL", "^", "AO L"}, {" :", "ABLE", "", "EY B AX L"}, {"", "ABLE", "", "AX B AX L"},
    {"", "ANG", "+", "EY N JH"}, {"", "A", "", "AE"},
    // B
    {" ", "BE", "^#", "B IH"}, {"", "BEING", "", "B IY IH NG"}, {" ", "BOTH", " ", "B OW TH"},
    {" ", "BUS", "#", "B IH Z"}, {"", "BUIL", "", "B IH L"}, {"", "B", "", "B"},
    // C
    {" ", "CH", "^", "K"}, {"^E", "CH", "", "K"}, {"", "CH", "", "CH"}, {" S", "CI", "#", "S AY"},
    {"", "CI", "A", "SH"}, {"", "CI", "O", "SH"}, {"", "CI", "EN", "SH"}, {"", "C", "+", "S"},
    {"", "CK", "", "K"}, {"", "COM", "%", "K AH M"}, {"", "C", "", "K"},
    // D
    {"#:", "DED", " ", "D IH D"}, {".E", "D", " ", "D"}, {"#^:E", "D", " ", "T"}, {" ", "DE", "^#", "D IH"},
    {" ", "DO", " ", "D UW"}, {" ", "DOES", "", "D AH Z"}, {" ", "DOING", "", "D UW IH NG"},
    {" ", "DOW", "", "D AW"}, {"", "DU", "A", "JH UW"}, {"", "D", "", "D"},
    // E
    {"#:", "E", " ", ""}, {"'^:", "E", " ", ""}, {" :", "E", " ", "IY"}, {"#", "ED", " ", "D"},
    {"#:", "E", "D ", ""}, {"", "EV", "ER", "EH V"}, {"", "E", "^%", "IY"}, {"", "ERI", "#", "IY R IY"},
    {"", "ERI", "", "EH R IH"}, {"#:", "ER", "#", "ER"}, {"", "ER", "#", "EH R"}, {"", "ER", "", "ER"},
    {" ", "EVEN", "", "IY V EH N"}, {"#:", "E", "W", ""}, {"@", "EW", "", "UW"}, {"", "EW", "", "Y UW"},
    {"", "E", "O", "IY"}, {"#:&", "ES", " ", "IH Z"}, {"#:", "E", "S ", ""}, {"#:", "ELY", " ", "L IY"},
    {"#:", "EMENT", "", "M EH N T"}, {"", "EFUL", "", "F UH L"}, {"", "EE", "", "IY"}, {"", "EARN", "", "ER N"},
    {" ", "EAR", "^", "ER"}, {"", "EAD", "", "EH D"}, {"#:", "EA", " ", "IY AX"}, {"", "EA", "SU", "EH"},
    {"", "EA", "", "IY"}, {"", "EIGH", "", "EY"}, {"", "EI", "", "IY"}, {" ", "EYE", "", "AY"},
    {"", "EY", "", "IY"}, {"", "EU", "", "Y UW"}, {"", "E", "", "EH"},
    // F
    {"", "FUL", "", "F UH L"}, {"", "F", "", "F"},
    // G
    {"", "GIV", "", "G IH V"}, {" ", "G", "I^", "G"}, {"", "GE", "T", "G EH"}, {"SU", "GGES", "", "G JH EH S"},
    {"", "GG", "", "G"}, {" B#", "G", "", "G"}, {"", "G", "+", "JH"}, {"", "GREAT", "", "G R EY T"},
    {"#", "GH", "", ""}, {"", "G", "", "G"},
    // H
    {" ", "HAV", "", "HH AE V"}, {" ", "HERE", "", "HH IY R"}, {" ", "HOUR", "", "AW ER"},
    {"", "HOW", "", "HH AW"}, {"", "H", "#", "HH"}, {"", "H", "", ""},
    // I
    {" ", "IN", "", "IH N"}, {" ", "I", " ", "AY"}, {"", "IN", "D", "AY N"}, {"", "IER", "", "IY ER"},
    {"#:R", "IED", " ", "IY D"}, {"", "IED", " ", "AY D"}, {"", "IEN", "", "IY EH N"}, {"", "IE", "T", "AY EH"},
    {" :", "I", "%", "AY"}, {"", "I", "%", "IY"}, {"", "IE", "", "IY"}, {"", "I", "^+:#", "IH"},
    {"", "IR", "#", "AY R"}, {"", "IZ", "%", "AY Z"}, {"", "IS", "%", "AY Z"}, {"", "I", "D%", "AY"},
    {"+^", "I", "^+", "IH"}, {"", "I", "T%", "AY"}, {"#^:", "I", "^+", "IH"}, {"", "I", "^+", "AY"},
    {"", "IR", "", "ER"}, {"", "IGH", "", "AY"}, {"", "ILD", "", "AY L D"}, {"", "IGN", " ", "AY N"},
    {"", "IGN", "^", "AY N"}, {"", "IGN", "%", "AY N"}, {"", "IQUE", "", "IY K"}, {"", "I", "", "IH"},
    // J K
    {"", "J", "", "JH"}, {" ", "K", "N", ""}, {"", "K", "", "K"},
    // L
    {"", "LO", "C#", "L OW"}, {"L", "L", "", ""}, {"#^:", "L", "%", "AX L"}, {"", "LEAD", "", "L IY D"},
    {"", "L", "", "L"},
    // M
    {"", "MOV", "", "M UW V"}, {"", "M", "", "M"},
    // N
    {"E", "NG", "+", "N JH"}, {"", "NG", "R", "NG G"}, {"", "NG", "#", "NG G"}, {"", "NGL", "%", "NG G AX L"},
    {"", "NG", "", "NG"}, {"", "NK", "", "NG K"}, {" ", "NOW", " ", "N AW"}, {"", "N", "", "N"},
    // O
    {"", "OF", " ", "AX V"}, {"", "OROUGH", "", "ER OW"}, {"#:", "OR", " ", "ER"}, {"#:", "ORS", " ", "ER Z"},
    {"", "OR", "", "AO R"}, {" ", "ONE", "", "W AH N"}, {"", "OW", "", "OW"}, {" ", "OVER", "", "OW V ER"},
    {"", "OV", "", "AH V"}, {"", "O", "^%", "OW"}, {"", "O", "^EN", "OW"}, {"", "O", "^I#", "OW"},
    {"", "OL", "D", "OW L"}, {"", "OUGHT", "", "AO T"}, {"", "OUGH", "", "AH F"}, {" ", "OU", "", "AW"},
    {"H", "OU", "S#", "AW"}, {"", "OUS", "", "AX S"}, {"", "OUR", "", "AO R"}, {"", "OULD", "", "UH D"},
    {"^", "OU", "^L", "AH"}, {"", "OUP", "", "UW P"}, {"", "OU", "", "AW"}, {"", "OY", "", "OY"},
    {"", "OING", "", "OW IH NG"}, {"", "OI", "", "OY"}, {"", "OOR", "", "AO R"}, {"", "OOK", "", "UH K"},
    {"", "OOD", "", "UH D"}, {"", "OO", "", "UW"}, {"", "O", "E", "OW"}, {"", "O", " ", "OW"},
    {"", "OA", "", "OW"}, {" ", "ONLY", "", "OW N L IY"}, {" ", "ONCE", "", "W AH N S"},
    {"", "ON'T", "", "OW N T"}, {"C", "O", "N", "AA"}, {"", "O", "NG", "AO"}, {" :^", "O", "N", "AH"},
    {"I", "ON", "", "AX N"}, {"#:", "ON", " ", "AX N"}, {"#^", "ON", "", "AX N"}, {"", "O", "ST ", "OW"},
    {"", "OF", "^", "AO F"}, {"", "OTHER", "", "AH DH ER"}, {"", "OSS", " ", "AO S"}, {"#^:", "OM", "", "AH M"},
    {"", "O", "", "AA"},
    // P Q R
    {"", "PH", "", "F"}, {"", "PEOP", "", "P IY P"}, {"", "POW", "", "P AW"}, {"", "PUT", " ", "P UH T"},
    {"", "P", "", "P"}, {"", "QUAR", "", "K W AO R"}, {"", "QU", "", "K W"}, {"", "Q", "", "K"},
    {" ", "RE", "^#", "R IY"}, {"", "R", "", "R"},
    // S
    {"", "SH", "", "SH"}, {"#", "SION", "", "ZH AX N"}, {"", "SOME", "", "S AH M"}, {"#", "SUR", "#", "ZH ER"},
    {"", "SUR", "#", "SH ER"}, {"#", "SU", "#", "ZH UW"}, {"#", "SSU", "#", "SH UW"}, {"#", "SED", " ", "Z D"},
    {"#", "S", "#", "Z"}, {"", "SAID", "", "S EH D"}, {"^", "SION", "", "SH AX N"}, {"", "S", "S", ""},
    {".", "S", " ", "Z"}, {"#:.E", "S", " ", "Z"}, {"#^:##", "S", " ", "Z"}, {"#^:#", "S", " ", "S"},
    {"U", "S", " ", "S"}, {" :#", "S", " ", "Z"}, {" ", "SCH", "", "S K"}, {"", "S", "C+", ""},
    {"#", "SM", "", "Z M"}, {"#", "SN", "'", "Z AX N"}, {"", "S", "", "S"},
    // T
    {" ", "THE", " ", "DH AX"}, {"", "TO", " ", "T UW"}, {"", "THAT", " ", "DH AE T"}, {" ", "THIS", " ", "DH IH S"},
    {" ", "THEY", "", "DH EY"}, {" ", "THERE", "", "DH EH R"}, {"", "THER", "", "DH ER"},
    {"", "THEIR", "", "DH EH R"}, {" ", "THAN", " ", "DH AE N"}, {" ", "THEM", " ", "DH EH M"},
    {"", "THESE", " ", "DH IY Z"}, {" ", "THEN", "", "DH EH N"}, {"", "THROUGH", "", "TH R UW"},
    {"", "THOSE", "", "DH OW Z"}, {"", "THOUGH", " ", "DH OW"}, {" ", "THUS", "", "DH AH S"},
    {"", "TH", "", "TH"}, {"#:", "TED", " ", "T IH D"}, {"S", "TI", "#N", "CH"}, {"", "TI", "O", "SH"},
    {"", "TI", "A", "SH"}, {"", "TIEN", "", "SH AX N"}, {"", "TUR", "#", "CH ER"}, {"", "TU", "A", "CH UW"},
    {" ", "TWO", "", "T UW"}, {"", "T", "", "T"},
    // U
    {" ", "UN", "I", "Y UW N"}, {" ", "UN", "", "AH N"}, {" ", "UPON", "", "AX P AO N"}, {"@", "UR", "#", "UH R"},
    {"", "UR", "#", "Y UH R"}, {"", "UR", "*", "ER"}, {"", "U", "^ ", "AH"}, {"", "U", "^^", "AH"},
    {"", "UY", "", "AY"}, {" G", "U", "#", ""}, {"G", "U", "%", ""}, {"G", "U", "#", "W"},
    {"#N", "U", "", "Y UW"}, {"@", "U", "", "UW"}, {"", "U", "", "Y UW"},
    // V W
    {"", "VIEW", "", "V Y UW"}, {"", "V", "", "V"}, {" ", "WERE", "", "W ER"}, {"", "WA", "S", "W AA"},
    {"", "WA", "T", "W AA"}, {"", "WHERE", "", "W EH R"}, {"", "WHAT", "", "W AA T"}, {"", "WHOL", "", "HH OW L"},
    {"", "WHO", "", "HH UW"}, {"", "WH", "", "W"}, {"", "WAR", "", "W AO R"}, {"", "WOR", "^", "W ER"},
    {"", "WR", "", "R"}, {"", "W", "", "W"},
    // X Y Z
    {"", "X", "", "K S"}, {"", "YOUNG", "", "Y AH NG"}, {" ", "YOU", "", "Y UW"}, {" ", "YES", "", "Y EH S"},
    {" ", "Y", "", "Y"}, {"#^:", "Y", " ", "IY"}, {"#^:", "Y", "I", "IY"}, {" :", "Y", " ", "AY"},
    {" :", "Y", "#", "AY"}, {" :", "Y", "^+:#", "IH"}, {" :", "Y", "^#", "AY"}, {"", "Y", "", "IH"},
    {"", "Z", "", "Z"},
};

bool is_vowel_letter(char letter) {
  return letter == 'A' || letter == 'E' || letter == 'I' || letter == 'O' || letter == 'U' || letter == 'Y';
}

bool is_consonant_letter(char letter) {
  return std::isalpha(static_cast<unsigned char>(letter)) != 0 && !is_vowel_letter(letter);
}

bool is_voiced_letter(char letter) {
  return std::strchr("BDVGJLMNRWZ", letter) != nullptr;
}

bool is_front_vowel_letter(char letter) {
  return letter == 'E' || letter == 'I' || letter == 'Y';
}

// Matches a right context starting at `position` (moving forward).
bool match_right(const std::string& word, int position, const char* pattern) {
  int index = position;
  int length = static_cast<int>(word.size());
  for (const char* cursor = pattern; *cursor != '\0'; ++cursor) {
    char symbol = *cursor;
    char letter = index < length ? word[static_cast<std::size_t>(index)] : ' ';
    switch (symbol) {
      case '#': {
        if (!is_vowel_letter(letter)) {
          return false;
        }
        while (index < length && is_vowel_letter(word[static_cast<std::size_t>(index)])) {
          ++index;
        }
        break;
      }
      case '*': {
        if (!is_consonant_letter(letter)) {
          return false;
        }
        while (index < length && is_consonant_letter(word[static_cast<std::size_t>(index)])) {
          ++index;
        }
        break;
      }
      case ':': {
        while (index < length && is_consonant_letter(word[static_cast<std::size_t>(index)])) {
          ++index;
        }
        break;
      }
      case '^':
        if (!is_consonant_letter(letter)) {
          return false;
        }
        ++index;
        break;
      case '.':
        if (!is_voiced_letter(letter)) {
          return false;
        }
        ++index;
        break;
      case '+':
        if (!is_front_vowel_letter(letter)) {
          return false;
        }
        ++index;
        break;
      case '$': {
        if (!is_consonant_letter(letter)) {
          return false;
        }
        char after = index + 1 < length ? word[static_cast<std::size_t>(index + 1)] : ' ';
        if (after != 'E' && after != 'I') {
          return false;
        }
        index += 2;
        break;
      }
      case '%': {
        const char* suffixes[6] = {"ER", "ES", "ED", "ING", "ELY", "E"};
        bool found = false;
        for (const char* suffix : suffixes) {
          int suffix_length = static_cast<int>(std::strlen(suffix));
          if (word.compare(static_cast<std::size_t>(std::min(index, length)), static_cast<std::size_t>(suffix_length),
                           suffix) == 0) {
            index += suffix_length;
            found = true;
            break;
          }
        }
        if (!found) {
          return false;
        }
        break;
      }
      case '&': {
        if (letter == 'S' || letter == 'C') {
          char after = index + 1 < length ? word[static_cast<std::size_t>(index + 1)] : ' ';
          index += after == 'H' ? 2 : 1;
        } else if (std::strchr("GZXJ", letter) != nullptr) {
          ++index;
        } else {
          return false;
        }
        break;
      }
      case '@': {
        if (letter == 'T' || letter == 'C' || letter == 'S') {
          char after = index + 1 < length ? word[static_cast<std::size_t>(index + 1)] : ' ';
          index += after == 'H' ? 2 : 1;
        } else if (std::strchr("RDLZNJ", letter) != nullptr) {
          ++index;
        } else {
          return false;
        }
        break;
      }
      case ' ':
        if (index < length && letter != ' ') {
          return false;
        }
        ++index;
        break;
      default:
        if (letter != symbol) {
          return false;
        }
        ++index;
        break;
    }
  }
  return true;
}

// Matches a left context ending just before `position` (moving backward).
bool match_left(const std::string& word, int position, const char* pattern) {
  int index = position - 1;
  int pattern_length = static_cast<int>(std::strlen(pattern));
  for (int cursor = pattern_length - 1; cursor >= 0; --cursor) {
    char symbol = pattern[cursor];
    char letter = index >= 0 ? word[static_cast<std::size_t>(index)] : ' ';
    switch (symbol) {
      case '#':
        if (!is_vowel_letter(letter)) {
          return false;
        }
        while (index >= 0 && is_vowel_letter(word[static_cast<std::size_t>(index)])) {
          --index;
        }
        break;
      case '*':
        if (!is_consonant_letter(letter)) {
          return false;
        }
        while (index >= 0 && is_consonant_letter(word[static_cast<std::size_t>(index)])) {
          --index;
        }
        break;
      case ':': {
        // Zero or more consonants, but when one consonant ('^') is required
        // next, leave it one: greedy here, "#^:" never matched "lazy" and
        // every final -y after a consonant came out /I/ instead of /i/.
        int run = 0;
        while (index - run >= 0 && is_consonant_letter(word[static_cast<std::size_t>(index - run)])) {
          ++run;
        }
        if (cursor > 0 && pattern[cursor - 1] == '^' && run > 0) {
          --run;
        }
        index -= run;
        break;
      }
      case '^':
        if (!is_consonant_letter(letter)) {
          return false;
        }
        --index;
        break;
      case '.':
        if (!is_voiced_letter(letter)) {
          return false;
        }
        --index;
        break;
      case '+':
        if (!is_front_vowel_letter(letter)) {
          return false;
        }
        --index;
        break;
      case '&':
        if (letter == 'H' && index >= 1 &&
            (word[static_cast<std::size_t>(index - 1)] == 'S' || word[static_cast<std::size_t>(index - 1)] == 'C')) {
          index -= 2;
        } else if (std::strchr("SCGZXJ", letter) != nullptr) {
          --index;
        } else {
          return false;
        }
        break;
      case '@':
        if (letter == 'H' && index >= 1 && std::strchr("TCS", word[static_cast<std::size_t>(index - 1)]) != nullptr) {
          index -= 2;
        } else if (std::strchr("TSRDLZNJ", letter) != nullptr) {
          --index;
        } else {
          return false;
        }
        break;
      case ' ':
        if (index >= 0 && letter != ' ') {
          return false;
        }
        --index;
        break;
      default:
        if (letter != symbol) {
          return false;
        }
        --index;
        break;
    }
  }
  return true;
}

std::vector<std::string> split_phones(const std::string& text) {
  std::vector<std::string> phones;
  std::istringstream stream(text);
  std::string phone;
  while (stream >> phone) {
    phones.push_back(phone);
  }
  return phones;
}

std::vector<std::string> letter_to_sound(const std::string& lower) {
  std::string word;
  for (char letter : lower) {
    if (std::isalpha(static_cast<unsigned char>(letter)) != 0 || letter == '\'') {
      word.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letter))));
    }
  }
  std::vector<std::string> phones;
  int position = 0;
  int length = static_cast<int>(word.size());
  while (position < length) {
    bool matched = false;
    for (const Rule& rule : kRules) {
      int match_length = static_cast<int>(std::strlen(rule.match));
      if (rule.match[0] != word[static_cast<std::size_t>(position)]) {
        continue;
      }
      if (word.compare(static_cast<std::size_t>(position), static_cast<std::size_t>(match_length), rule.match) != 0) {
        continue;
      }
      if (!match_left(word, position, rule.left)) {
        continue;
      }
      if (!match_right(word, position + match_length, rule.right)) {
        continue;
      }
      std::vector<std::string> produced = split_phones(rule.output);
      phones.insert(phones.end(), produced.begin(), produced.end());
      position += match_length;
      matched = true;
      break;
    }
    if (!matched) {
      ++position;
    }
  }
  // Doubled letters spell one consonant ("happy" is HH AE P IY, not P P).
  std::vector<std::string> single;
  for (const std::string& phone : phones) {
    if (!single.empty() && single.back() == phone && !is_vowel_phone(phone)) {
      continue;
    }
    single.push_back(phone);
  }
  return single;
}

// Stress for rule-derived words: primary on the first full vowel, except
// after common unstressed prefixes, or before stress-attracting suffixes.
void assign_stress(const std::string& lower, std::vector<Phone>& phones) {
  std::vector<std::size_t> vowels;
  for (std::size_t index = 0; index < phones.size(); ++index) {
    if (is_vowel_phone(phones[index].name)) {
      vowels.push_back(index);
      phones[index].stress = 0;
    }
  }
  if (vowels.empty()) {
    return;
  }
  std::size_t chosen = 0;
  const char* prefixes[] = {"be", "de", "re", "con", "com", "ex", "pre", "pro", "a", "un", "dis", "mis", "in", "en"};
  for (const char* prefix : prefixes) {
    std::size_t prefix_length = std::strlen(prefix);
    if (vowels.size() >= 2 && lower.size() > prefix_length + 2 && lower.compare(0, prefix_length, prefix) == 0) {
      chosen = 1;
      break;
    }
  }
  const char* before_last[] = {"tion", "sion", "ic", "ical", "ity", "ial", "ian", "ious"};
  for (const char* suffix : before_last) {
    std::size_t suffix_length = std::strlen(suffix);
    if (lower.size() > suffix_length + 1 && lower.compare(lower.size() - suffix_length, suffix_length, suffix) == 0 &&
        vowels.size() >= 2) {
      chosen = vowels.size() - 2;
      if (std::strcmp(suffix, "ity") == 0 && vowels.size() >= 3) {
        chosen = vowels.size() - 3;
      }
      break;
    }
  }
  // A schwa is never the stressed vowel if a full vowel exists.
  if (phones[vowels[chosen]].name == "AX") {
    for (std::size_t index = 0; index < vowels.size(); ++index) {
      if (phones[vowels[index]].name != "AX") {
        chosen = index;
        break;
      }
    }
  }
  phones[vowels[chosen]].stress = 1;
  for (std::size_t index = 0; index < vowels.size(); ++index) {
    if (index != chosen && index + 2 <= chosen) {
      phones[vowels[index]].stress = 2;
    }
  }
}

class Lexicon {
 public:
  explicit Lexicon(const char* source) {
    std::istringstream stream(source);
    read(stream);
  }

  // A pronouncing dictionary in CMUdict form: "word PH ON1 ES", alternates
  // as "word(2)" (skipped: the first is the common one), "#" comments.
  explicit Lexicon(std::istream& stream) { read(stream); }

  bool empty() const { return entries_.empty(); }

  bool find(const std::string& word, std::vector<Phone>& phones) const {
    std::map<std::string, std::vector<Phone>>::const_iterator found = entries_.find(word);
    if (found == entries_.end()) {
      return false;
    }
    phones = found->second;
    return true;
  }

 private:
  void read(std::istream& stream) {
    std::string line;
    while (std::getline(stream, line)) {
      std::string::size_type comment = line.find('#');
      if (comment != std::string::npos) {
        line.erase(comment);
      }
      std::istringstream fields(line);
      std::string word;
      if (!(fields >> word) || word.find('(') != std::string::npos) {
        continue;
      }
      std::vector<Phone> phones;
      std::string token;
      while (fields >> token) {
        Phone phone;
        char last = token.back();
        if (last >= '0' && last <= '2') {
          phone.stress = last - '0';
          token.pop_back();
        }
        phone.name = token;
        phones.push_back(phone);
      }
      entries_[word] = phones;
    }
  }

  std::map<std::string, std::vector<Phone>> entries_;
};

const Lexicon& lexicon() {
  static const Lexicon instance(kLexicon);
  return instance;
}

// Weak forms: in running speech a function word's vowel reduces to schwa
// (or syllabic r). Aligned against the bdl labels, our unstressed "to"
// vowel was a schwa there 19 times in 27, and unstressed AE and EH ("and",
// "them") 12 times each, while we spoke the full vowels.
constexpr const char* kWeakForms = R"(a AX0
an AX0 N
the DH AX0
of AX0 V
to T AX0
and AX0 N D
at AX0 T
as AX0 Z
can K AX0 N
for F ER0
from F R AX0 M
was W AX0 Z
were W ER0
them DH AX0 M
than DH AX0 N
but B AX0 T
some S AX0 M
had HH AX0 D
has HH AX0 Z
have HH AX0 V
would W AX0 D
could K AX0 D
should SH AX0 D
her HH ER0
your Y ER0
)";

// The CMU Pronouncing Dictionary, when present beside the voice (ref/,
// fetched separately): our own lexicon is consulted first, the letter rules
// last. The rules alone mispronounced much of ARCTIC (42% of its word
// tokens went through them: "Philip" /faIlIp/, "delicate" stressed late).
const Lexicon& dictionary() {
  static std::ifstream file("ref/cmudict/cmudict.dict");
  static const Lexicon instance(file);
  return instance;
}

const Lexicon& weak_forms() {
  static const Lexicon instance(kWeakForms);
  return instance;
}

}  // namespace

bool is_vowel_phone(const std::string& name) {
  const char* vowels[] = {"AA", "AE", "AH", "AO", "AW", "AX", "AY", "EH", "ER", "EY", "IH", "IY", "OW", "OY", "UH", "UW"};
  for (const char* vowel : vowels) {
    if (name == vowel) {
      return true;
    }
  }
  return false;
}

bool weak_form(const std::string& word, std::vector<Phone>& phones) {
  return weak_forms().find(word, phones);
}

bool is_function_word(const std::string& word) {
  for (const char* entry : kFunctionWords) {
    if (word == entry) {
      return true;
    }
  }
  return false;
}

std::vector<Phone> pronounce(const std::string& word) {
  std::vector<Phone> phones;
  if (lexicon().find(word, phones)) {
    return phones;
  }
  // Plural and past forms of lexicon words.
  if (word.size() > 2 && word.back() == 's' && lexicon().find(word.substr(0, word.size() - 1), phones)) {
    std::string last = phones.back().name;
    Phone ending;
    bool voiceless = last == "P" || last == "T" || last == "K" || last == "F" || last == "TH";
    bool sibilant = last == "S" || last == "Z" || last == "SH" || last == "ZH" || last == "CH" || last == "JH";
    if (sibilant) {
      Phone vowel;
      vowel.name = "AX";
      vowel.stress = 0;
      phones.push_back(vowel);
    }
    ending.name = voiceless ? "S" : "Z";
    phones.push_back(ending);
    return phones;
  }
  if (dictionary().find(word, phones)) {
    return phones;
  }
  std::vector<std::string> names = letter_to_sound(word);
  for (const std::string& name : names) {
    Phone phone;
    phone.name = name;
    phones.push_back(phone);
  }
  assign_stress(word, phones);
  // An unstressed /I/ inside a rule-derived word is a schwa ("Philip",
  // "delicate", "capital": bdl's labels have ax where ours had ih, 24
  // times). The last vowel keeps /I/ ("needed", "chances" end in -Id, -Iz).
  std::size_t last_vowel = phones.size();
  for (std::size_t index = 0; index < phones.size(); ++index) {
    if (is_vowel_phone(phones[index].name)) {
      last_vowel = index;
    }
  }
  for (std::size_t index = 0; index < phones.size(); ++index) {
    if (phones[index].name == "IH" && phones[index].stress == 0 && index != last_vowel) {
      phones[index].name = "AX";
    }
  }
  return phones;
}

}  // namespace voice
