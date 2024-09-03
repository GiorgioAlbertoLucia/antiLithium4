root -b <<-EOF
auto inFile = TFile::Open("AO2D.root");
inFile->ls();
const int nKeys = inFile->GetNkeys();
TIter iterNext(inFile->GetListOfKeys());

TKey *key;
//while ((key = (TKey*)iterNext())) 
//{
//    std::cout << key->GetName() << ", " << key->GetClassName() << std::endl;
//}
EOF