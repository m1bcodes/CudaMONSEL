#ifndef _REFERENCE_CUH_
#define _REFERENCE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Reference
{
   class Reference
   {
   protected:
      Reference::Reference();

   public:
      virtual StringT getShortForm() const = 0;
      virtual StringT getLongForm() const = 0;
   };

   class Author
   {
   public:
      Author(StringT first, StringT last, StringT affiliation);
      Author(StringT first, StringT last);
      Author(const Author&);
      StringT getAuthor();
      StringT getAffiliation();
      StringT getFirst();
      StringT getLast();

   private:
      StringT mFirst;
      StringT mLast;
      StringT mAffiliation;
   };

   typedef std::vector<Author> AuthorList;

   static StringT toString(AuthorList authors);
   static StringT toAbbrev(AuthorList authors);

   class Journal
   {
   public:
      Journal(StringT name, StringT abbrev, StringT publisher);
      Journal(const Journal&);
      StringT getName();
      StringT getPublisher();

   private:
      StringT mName;
      StringT mPublisher;
      StringT mAbbreviation;
   };

   class JournalArticle : public Reference
   {
   public:
      JournalArticle(StringT title, Journal journal, StringT vol, StringT pages, int year, const Author authors[], int len);
      JournalArticle(const Journal& journal, StringT vol, StringT pages, int year, const Author authors[], int len);
      StringT getShortForm() const override;
      StringT getLongForm() const override;

   private:
      StringT mTitle;
      const Journal& mJournal;
      StringT mVolume;
      StringT mPages;
      int mYear;
      AuthorList mAuthors;
   };

   class Book : public Reference
   {
   private:
      StringT mTitle;
      StringT mPublisher;
      int mYear;
      AuthorList mAuthors;

   public:
      Book(StringT title, StringT publisher, int year, Author authors[], int len);
      StringT getShortForm() const override;
      StringT getLongForm() const override;
   };

   class Program : public Reference
   {
   private:
      StringT mTitle;
      StringT mVersion;
      AuthorList mAuthors;

   public:
      Program(StringT title, StringT version, Author authors[], int len);

      StringT getShortForm() const;
      StringT getLongForm() const;
   };

   class WebSite : public Reference
   {
   private:
      const StringT mUrl;
      const StringT mTitle;
      const StringT mDate;
      const AuthorList mAuthors;

   public:
      WebSite(StringT url, StringT title, StringT date, const Author authors[], int len);

      StringT getShortForm() const override;
      StringT getLongForm() const override;
   };

   class CrudeReference : public Reference
   {
   public:
      CrudeReference(StringT ref);

      StringT getShortForm() const override;
      StringT getLongForm() const override;

   private:
      StringT mReference;
   };

   extern const StringT ONERA;
   extern const StringT NIST;
   extern const StringT NBS;
   extern const StringT LehighU;
   extern const StringT Eindhoven;

   // Some prolific authors
   extern const Author DNewbury;
   extern const Author KHeinrich;
   extern const Author JPouchou;
   extern const Author FPichoir;
   extern const Author RCastaing;
   extern const Author RMyklebust;
   extern const Author DBright;
   extern const Author CFiori;
   extern const Author JArmstrong;
   extern const Author JSmall;
   extern const Author JGoldstein;
   extern const Author DWilliams;
   extern const Author GBastin;
   extern const Author HHeijligers;
   extern const Author RPackwood;
   extern const Author CLyman;
   extern const Author ELifshin;
   extern const Author PEchlin;
   extern const Author LSawyer;
   extern const Author DJoy;
   extern const Author JMichael;
   extern const Author PDuncumb;
   extern const Author PStatham;
   extern const Author SReed;
   extern const Author JPhilibert;
   extern const Author HYakowitz;
   extern const Author JCriss;
   extern const Author LBirks;
   extern const Author TMcKinley;
   extern const Author DWittry;
   extern const Author GLove;
   extern const Author VScott;
   extern const Author JDijkstra;
   extern const Author Savitzky;
   extern const Author Golay;
   extern const Author RFEdgerton;
   extern const Author Cullen;
   extern const Author CPowell;
   extern const Author FSalvat;
   extern const Author AJablonski;
   extern const Author BHenke;
   extern const Author EGullikson;
   extern const Author JDavis;
   extern const Author CSwytThomas;
   extern const Author VanGrieken;
   extern const Author Markowicz;
   extern const Author Oberndorff;
   extern const Author Czyzewski;
   extern const Author MacCallum;

   // Commonly referenced journals
   extern const Journal MandM;
   extern const Journal Scanning;
   extern const Journal Ultramicroscopy;
   extern const Journal SIA;
   extern const Journal JPhysD;
   extern const Journal XRaySpec;
   extern const Journal AnalChem;
   extern const Journal JApplPhys;
   extern const Journal AtDatNucData;
   extern const Journal PhysRevA;
   extern const Journal ApplPhysLett;
}

#endif