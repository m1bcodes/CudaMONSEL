#ifndef _REFERENCE_CUH_
#define _REFERENCE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "Amphibian\vector.cuh"

#include <cuda_runtime.h>

namespace Reference
{
   class Reference
   {
   protected:
      __host__ __device__ Reference::Reference();

   public:
      __host__ __device__ virtual StringT getShortForm() const = 0;
      __host__ __device__ virtual StringT getLongForm() const = 0;
   };

   class Author
   {
   public:
      __host__ __device__ Author(StringT first, StringT last, StringT affiliation);
      __host__ __device__ Author(StringT first, StringT last);
      __host__ __device__ Author(const Author&);
      __host__ __device__ StringT getAuthor() const;
      __host__ __device__ StringT getAffiliation() const;
      __host__ __device__ StringT getFirst() const;
      __host__ __device__ StringT getLast() const;

   private:
      StringT mFirst;
      StringT mLast;
      StringT mAffiliation;
   };

   typedef amp::vector<const Author*> AuthorList;

   class Journal
   {
   public:
      __host__ __device__ Journal(StringT name, StringT abbrev, StringT publisher);
      __host__ __device__ Journal(const Journal&);
      __host__ __device__ StringT getName();
      __host__ __device__ StringT getPublisher();

   private:
      StringT mName;
      StringT mPublisher;
      StringT mAbbreviation;
   };

   class JournalArticle : public Reference
   {
   public:
      __host__ __device__ JournalArticle(StringT title, const Journal& journal, StringT vol, StringT pages, int year, const Author* authors[], int len);
      __host__ __device__ JournalArticle(const Journal& journal, StringT vol, StringT pages, int year, const Author* authors[], int len);
      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

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
   public:
      __host__ __device__ Book(StringT title, StringT publisher, int year, const Author* authors[], int len);
      __host__ __device__ StringT GetTitle() const;
      __host__ __device__ StringT GetPublisher() const;
      __host__ __device__ int GetYear() const;
      __host__ __device__ const AuthorList& GetAuthors() const;

      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

   private:
      StringT mTitle;
      StringT mPublisher;
      int mYear;
      AuthorList mAuthors;
   };

   class Program : public Reference
   {
   public:
      __host__ __device__ Program(StringT title, StringT version, const Author* authors[], int len);

      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

   private:
      StringT mTitle;
      StringT mVersion;
      AuthorList mAuthors;
   };

   class BookChapter : public Reference
   {
   public:
      __host__ __device__ BookChapter(const Book& book, StringT pages, const Author* authors[], int len);
      __host__ __device__ BookChapter(const Book& book, const Author* authors[], int len);

      __host__ __device__ StringT getPages() const;
      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

   private:
      Book mBook;
      StringT mPages;
      AuthorList mAuthors;
   };

   class WebSite : public Reference
   {
   public:
      __host__ __device__ WebSite(StringT url, StringT title, StringT date, const Author* authors[], int len);

      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

   private:
      const StringT mUrl;
      const StringT mTitle;
      const StringT mDate;
      const AuthorList mAuthors;
   };

   class CrudeReference : public Reference
   {
   public:
      __host__ __device__ CrudeReference(StringT ref);

      __host__ __device__ StringT getShortForm() const override;
      __host__ __device__ StringT getLongForm() const override;

      __host__ __device__ const StringT &getReference() const;
      __host__ __device__ void setReference(const StringT& newref);

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
   extern const Author Bote;

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

   // book
   extern const Book ElectronProbeQuant;
   extern const Book GoldsteinBook;
   extern const Book QuantitativeElectronProbeMicroanalysisBook;
   extern const Book ElectronBeamXRayMicroanalysis;
   extern const Book EnergyDispersiveXRaySpectrometery;
   extern const Book CharacterizationOfParticles;
   extern const Book FrameBook;
   extern const Book ElectronMicroprobe;
   extern const Book MicrobeamAnalysis;
   extern const Book HandbookOfXRaySpectrometry;

   // default
   extern const CrudeReference NullReference;
   extern __device__ CrudeReference *d_NullReference;

   extern __global__ void initCuda(char *d_data);

}

#endif