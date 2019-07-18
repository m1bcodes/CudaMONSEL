#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"
#include "CudaUtil.h"

namespace Reference
{
   __host__ __device__ Reference::Reference()
   {
   }

   __host__ __device__ Author::Author(StringT first, StringT last, StringT affiliation) : mFirst(first), mLast(last), mAffiliation(affiliation)
   {
   }

   __host__ __device__ Author::Author(StringT first, StringT last) : mFirst(first), mLast(last)
   {
   }

   __host__ __device__ Author::Author(const Author& other) : mFirst(other.mFirst), mLast(other.mLast), mAffiliation(other.mAffiliation)
   {
   }

   __host__ __device__ StringT Author::getAuthor() const
   {
      return mLast + " " + mFirst.substr(0, 1);
   }

   __host__ __device__ StringT Author::getAffiliation() const
   {
      return mAffiliation;
   }

   __host__ __device__ StringT Author::getFirst() const
   {
      return mFirst;
   }

   __host__ __device__ StringT Author::getLast() const
   {
      return mLast;
   }

   __host__ __device__ StringT toString(const AuthorList& authors)
   {
      StringT ret;
      for (int i = 0; i < authors.size(); ++i) {
         if (i != 0)
            ret += (i == authors.size() - 1 ? " & " : ", ");
         ret += (authors.at(i)->getAuthor());
      }
      return ret;
   }

   __host__ __device__ StringT toAbbrev(const AuthorList& authors)
   {
      if (authors.size() > 1)
         return authors.at(0)->getLast() + " et al";
      else if (authors.size() > 0)
         return authors.at(0)->getLast();
      else
         return "";
   }

   // Some prolific establishments
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
#else
   const StringT ONERA("Office National d'Etudes et de Recherche Aerospatiales");
   const StringT NIST("National Institute of Standards & Technology");
   const StringT NBS("National Bureau of Standards");
   const StringT LehighU("Lehigh University");
   const StringT Eindhoven("University of Technology, Eindhoven");

   // Some prolific authors
   const Author DNewbury("Dale", "Newbury", NIST);
   const Author KHeinrich("Kurt", "Heinrich", NIST);
   const Author JPouchou("Jean-Louis", "Pouchou", ONERA);
   const Author FPichoir("Franï¿½iose", "Pichoir", ONERA);
   const Author RCastaing("R", "Castaing", "Universitï¿?de Paris-Sud");
   const Author RMyklebust("Robert", "Myklebust", NIST);
   const Author DBright("David", "Bright", NIST);
   const Author CFiori("Chuck", "Fiori", NBS);
   const Author JArmstrong("John", "Armstrong", "American University");
   const Author JSmall("John", "Small", NIST);
   const Author JGoldstein("Joseph", "Goldstein", LehighU);
   const Author DWilliams("Dave", "Williams", LehighU);
   const Author GBastin("G", "Bastin", Eindhoven);
   const Author HHeijligers("H", "Heijligers", Eindhoven);
   const Author RPackwood("Rod", "Packwood", "Metals Technology Laboratory");
   const Author CLyman("Charles", "Lyman", LehighU);
   const Author ELifshin("Eric", "Lifshin", "State University of New York at Albany");
   const Author PEchlin("Patrick", "Echlin", "Cambridge Analytical Microscopy, Ltd.");
   const Author LSawyer("Linda", "Sawyer", "Ticona, LLC");
   const Author DJoy("David", "Joy", "University of Tennessee");
   const Author JMichael("Joseph", "Michael", "Sandia National Laboratories");
   const Author PDuncumb("Peter", "Duncumb", "University of Cambridge");
   const Author PStatham("Peter", "Statham", "Oxford Instruments");
   const Author SReed("S. J.", "Reed", "");
   const Author JPhilibert("J.", "Philibert", "");
   const Author HYakowitz("Harvey", "Yakowitz", NBS);
   const Author JCriss("J.", "Criss", "");
   const Author LBirks("L. S.", "Birks", "");
   const Author TMcKinley("T. D.", "McKinley", "");
   const Author DWittry("D. B.", "Wittry", "");
   const Author GLove("G.", "Love", "");
   const Author VScott("V. D.", "Scott", "");
   const Author JDijkstra("J. M.", "Dijkstra", "");
   const Author Savitzky("A", "Savitzky");
   const Author Golay("M. J. E.", "Golay");
   const Author RFEdgerton("R. J.", "Edgerton");
   const Author Cullen("Dermott", "Cullen", "Lawrence Livermore National Laboratory");
   const Author CPowell("Cedric", "Powell", "N.I.S.T.");
   const Author FSalvat("Francesc", "Salvat", "Facultat de FÃ­sica (ECM), Universitat de Barcelona, Diagonal 647, 08028 Barcelona, Spain");
   const Author AJablonski("A", "Jablonksi");
   const Author BHenke("B.L.", "Henke");
   const Author EGullikson("E.M.", "Gullikson");
   const Author JDavis("J.C.", "Davis");
   const Author CSwytThomas("C.", "Swyt-Thomas", "N.I.H.");
   const Author VanGrieken("Rene", "Van Grieken", "");
   const Author Markowicz("Andrezej", "Markowicz", "");
   const Author Oberndorff("P. J. T. L.", "Oberndorff", "");
   const Author Czyzewski("", "Czyzewski", "");
   const Author MacCallum("", "MacCallum", "");
   const Author Bote("David", "Bote", "Facultat de Física (ECM), Universitat de Barcelona, Diagonal 647, 08028 Barcelona, Spain");
#endif

   __host__ __device__ Journal::Journal(StringT name, StringT abbrev, StringT publisher) : mName(name), mAbbreviation(abbrev), mPublisher(publisher)
   {
   }

   __host__ __device__ Journal::Journal(const Journal& other) : mName(other.mName), mAbbreviation(other.mAbbreviation), mPublisher(other.mPublisher)
   {
   }

   __host__ __device__ StringT Journal::getName()
   {
      return mName;
   }

   __host__ __device__ StringT Journal::getPublisher()
   {
      return mPublisher;
   }

   // Commonly referenced journals
   const Journal MandM("Microscopy and Microanalysis", "Microsc. Microanal. (USA)", "Cambridge University Press");
   const Journal Scanning("Scanning", "Scanning (USA)", "Foundation for the Advancement of Medicine & Science");
   const Journal Ultramicroscopy("Ultramicroscopy", "Ultramicroscopy (Netherlands)", "Elsevier Science");
   const Journal SIA("Surface and Interface Analysis", "Surf. Interface Anal. (UK)", "John Wiley & Sons Ltd.");
   const Journal JPhysD("Journal of Physics - D", "J. Phys. D.", "");
   const Journal XRaySpec("X-Ray Spectrometry", "X-Ray Spectrom.", "John Wiley & Sons Ltd.");
   const Journal AnalChem("Analytical Chemistry", "Anal Chem", "American Chemical Society");
   const Journal JApplPhys("Journal of Applied Physics", "J. Appl. Phys", "American Physical Society");
   const Journal AtDatNucData("Atomic Data and Nuclear Data Tables", "At. Dat. Nucl. Dat. Tables", "Academic Press");
   const Journal PhysRevA("Physical Review A", "Phys. Rev. A", "American Physical Society");
   const Journal ApplPhysLett("Applied Physics Letters", "Appl. Phys. Let.", "American Physical Society");

   __host__ __device__ JournalArticle::JournalArticle(StringT title, const Journal& journal, StringT vol, StringT pages, int year, const Author* authors[], int len) : Reference(), mAuthors(authors, authors + len), mJournal(journal), mTitle(title), mVolume(vol), mPages(pages), mYear(year)
   {
   }

   __host__ __device__ JournalArticle::JournalArticle(const Journal& journal, StringT vol, StringT pages, int year, const Author* authors[], int len) : mJournal(journal), mAuthors(authors, authors + len), mVolume(vol), mPages(pages), mYear(year)
   {
   }

   __host__ __device__ StringT JournalArticle::getShortForm() const
   {
      //return String.format("%s. %s %s p%s (%d)", Author.toAbbrev(mAuthors), mJournal.mAbbreviation, mVolume, mPages, mYear);
      char yearstr[8];
      amp::IToA(mYear, yearstr);
      return toAbbrev(mAuthors) + " " + mVolume + " " + mPages + " " + yearstr; //missing abbrev
   }

   __host__ __device__ StringT JournalArticle::getLongForm() const
   {
      //return Author::toString(mAuthors) + " " + mTitle + " " + mJournal.mAbbreviation + " " + mVolume + std::to_string(mPages) + " " + std::to_string(mYear);
      char yearstr[8];
      amp::IToA(mYear, yearstr);
      return toString(mAuthors) + " " + mTitle + " " + mVolume + mPages + " " + yearstr;
   }

   const Author* auLoveScott1978[] = { &GLove, &VScott };
   JournalArticle LoveScott1978(JPhysD, "11", "p 1369", 1978, auLoveScott1978, 2);

   const Author* alProza96[] = { &GBastin, &JDijkstra, &HHeijligers };
   JournalArticle Proza96(XRaySpec, "27", "p 3-10", 1998, alProza96, 3);

   const Author* alProza96Extended[] = { &GBastin, &Oberndorff, &JDijkstra, &HHeijligers };
   JournalArticle Proza96Extended(XRaySpec, "30", "p 382-387", 2001, alProza96Extended, 4);

   __host__ __device__ Book::Book(StringT title, StringT publisher, int year, const Author* authors[], int len) : mTitle(title), mPublisher(publisher), mYear(year), mAuthors(authors, authors + len)
   {
   }

   __host__ __device__ StringT Book::GetTitle() const
   {
      return mTitle;
   }

   __host__ __device__ StringT Book::GetPublisher() const
   {
      return mPublisher;
   }

   __host__ __device__ int Book::GetYear() const
   {
      return mYear;
   }

   __host__ __device__ const AuthorList& Book::GetAuthors() const
   {
      return mAuthors;
   }

   __host__ __device__ StringT Book::getShortForm() const
   {
      char yearstr[8];
      amp::IToA(mYear, yearstr);
      return toAbbrev(mAuthors) + mTitle + mPublisher + yearstr;
   }

   __host__ __device__ StringT Book::getLongForm() const
   {
      char yearstr[8];
      amp::IToA(mYear, yearstr);
      return toString(mAuthors) + mTitle + mPublisher + yearstr;
   }

   __host__ __device__ Program::Program(StringT title, StringT version, const Author* authors[], int len) : mTitle(title), mVersion(version), mAuthors(authors, authors + len)
   {
   }

   __host__ __device__ StringT Program::getShortForm() const
   {
      return toAbbrev(mAuthors) + mTitle + mVersion;
   }

   __host__ __device__ StringT Program::getLongForm() const
   {
      return toString(mAuthors) + mTitle + mVersion;
   }

   const Author* alENDLIB97_Relax[] = { &Cullen };
   const Program ENDLIB97_Relax("RELAX", "ENDLIB97", alENDLIB97_Relax, 1);

   const Author* alDTSA[] = { &CFiori, &CSwytThomas, &RMyklebust };
   const Program DTSA("DTSA", "3.0.1", alDTSA, 3);

   const Author* alElectronProbeQuant[] = { &KHeinrich, &DNewbury };
   const Book ElectronProbeQuant("Electron Probe Quantitation", "Plenum", 1991, alElectronProbeQuant, 2);

   const Author* alGoldsteinBook[] = { &JGoldstein, &CLyman, &DNewbury, &ELifshin, &PEchlin, &LSawyer, &DJoy, &JMichael };
   const Book GoldsteinBook("Scanning Electron Microscopy and X-Ray Microanalysis - 3rd edition", "Kluwer Academic/Plenum", 2003, alGoldsteinBook, 8);

   const Author* alQuantitativeElectronProbeMicroanalysisBook[] = { &KHeinrich };
   const Book QuantitativeElectronProbeMicroanalysisBook("Quantitative Electron Probe Microanalysis - NBS SP 298", "National Bureau of Standards", 1968, alQuantitativeElectronProbeMicroanalysisBook, 1);

   const Author* alElectronBeamXRayMicroanalysis[] = { &KHeinrich };
   const Book ElectronBeamXRayMicroanalysis("Electron Beam X-Ray Microanalysis", "Van Nostrand Reinhold Company", 1981, alElectronBeamXRayMicroanalysis, 1);

   const Author* alEnergyDispersiveXRaySpectrometery[] = { &KHeinrich, &DNewbury, &RMyklebust, &CFiori };
   const Book EnergyDispersiveXRaySpectrometery("Energy Dispersive X-Ray Spectrometery - NBS SP 604", "National Bureau of Standards", 1981, alEnergyDispersiveXRaySpectrometery, 4);

   const Author* alCharacterizationOfParticles[] = {&KHeinrich};
   const Book CharacterizationOfParticles("Characterization of Particles NBS SP 533", "National Bureau of Standards", 1978, alCharacterizationOfParticles, 1);

   const Author* alFrameBook[] = { &HYakowitz, &RMyklebust, &KHeinrich };
   const Book FrameBook("FRAME: An On-Line Correction Procedure for Quantitative Electron Probe Micronanalysis, NBS Tech Note 796", "National Bureau of Standards", 1973, alFrameBook, 3);

   const Author* alElectronMicroprobe[] = { &TMcKinley, &KHeinrich, &DWittry };
   const Book ElectronMicroprobe("The Electron Microprobe", "Wiley (New York)", 1966, alElectronMicroprobe, 3);

   const Author* alMicrobeamAnalysis[] = { &KHeinrich };
   const Book MicrobeamAnalysis("Microbeam Analysis", "San Francisco Press", 1982, alMicrobeamAnalysis, 1);

   const Author* alHandbookOfXRaySpectrometry[] { &VanGrieken, &Markowicz };
   const Book HandbookOfXRaySpectrometry("Handbook of X-Ray Spectrometry", "Marcel Dekker", 2002, alHandbookOfXRaySpectrometry, 2);

   const Author* alHenke1993[] { &BHenke, &EGullikson, &JDavis };
   const JournalArticle Henke1993("X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92", AtDatNucData, "54", "181-342", 1993, alHenke1993, 3);

   __host__ __device__ BookChapter::BookChapter(const Book& book, StringT pages, const Author* authors[], int len) : mBook(book), mPages(pages), mAuthors(authors, authors + len)
   {
   }

   __host__ __device__ BookChapter::BookChapter(const Book& book, const Author* authors[], int len) : mBook(book), mAuthors(authors, authors + len)
   {
   }

   __host__ __device__ StringT BookChapter::getPages() const
   {
      return mPages;
   }

   __host__ __device__ StringT BookChapter::getShortForm() const
   {
      char yearstr[8];
      amp::IToA(mBook.GetYear(), yearstr);
      return toAbbrev(mAuthors) + mBook.GetTitle() + toAbbrev(mBook.GetAuthors()) + mBook.GetPublisher() + yearstr;
   }

   __host__ __device__ StringT BookChapter::getLongForm() const
   {
      char yearstr[8];
      amp::IToA(mBook.GetYear(), yearstr);
      return toString(mAuthors) + mBook.GetTitle() + toString(mBook.GetAuthors()) + mBook.GetPublisher() + yearstr;
   }

   const Author* AuthorPAPinEPQ[] = { &JPouchou, &FPichoir };
   const BookChapter PAPinEPQ(ElectronProbeQuant, "p XXX-XXX", AuthorPAPinEPQ, 2);

   const Author* alJTA1982[] = { &JArmstrong };
   const BookChapter JTA1982(MicrobeamAnalysis, "175-180", alJTA1982, 1);

   //Date createDate(int yr, int month, int day)
   //{
   //   final Calendar cal = Calendar.getInstance(Locale.US);
   //   cal.set(yr, month, 24);
   //   return cal.getTime();
   //}

   __host__ __device__ WebSite::WebSite(StringT url, StringT title, StringT date, const Author* authors[], int len) : mUrl(url), mTitle(title), mDate(date), mAuthors(authors, authors + len)
   {
   }

   __host__ __device__ StringT WebSite::getShortForm() const
   {
      return mUrl;
   }

   __host__ __device__ StringT WebSite::getLongForm() const
   {
      return toString(mAuthors) + ". " + mTitle + "[" + mUrl + " on " + mDate + "]";
   }

   __host__ __device__ CrudeReference::CrudeReference(StringT ref) : mReference(ref)
   {
   }

   __host__ __device__ StringT CrudeReference::getShortForm() const
   {
      return mReference;
   }

   __host__ __device__ StringT CrudeReference::getLongForm() const
   {
      return mReference;
   }

   __host__ __device__ const StringT &CrudeReference::getReference() const
   {
      return mReference;
   }

   __host__ __device__ void CrudeReference::setReference(const StringT& newref)
   {
      mReference = newref;
   }

   const CrudeReference NullReference("-");
   __device__ CrudeReference *d_NullReference;

   __global__ void initCuda(char *d_data)
   {
      d_NullReference = new CrudeReference(d_data);
   }
}
