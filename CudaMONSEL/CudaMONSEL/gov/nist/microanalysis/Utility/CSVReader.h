#ifndef _CSVREADER_H_
#define _CSVREADER_H_

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class CSVRow
{
public:
   std::string const& operator[](std::size_t index) const;
   std::size_t size() const;
   void readNextRow(std::istream& str);

private:
   std::vector<std::string> m_data;
};

//int main()
//{
//   std::ifstream file("plop.csv");
//
//   CSVRow row;
//   while (file >> row)
//   {
//      std::cout << "4th Element(" << row[3] << ")\n";
//   }
//   file.close();
//}

class CSVIterator
{
public:
   //typedef std::input_iterator_tag iterator_category;
   //typedef CSVRow value_type;
   //typedef std::size_t difference_type;
   //typedef CSVRow* pointer;
   //typedef CSVRow& reference;

   static bool IsNaN(std::string a);

   CSVIterator(std::istream& str);
   CSVIterator();

   // Pre Increment
   CSVIterator& operator++();

   // Post increment
   CSVIterator operator++(int);

   CSVRow const& operator*() const;
   CSVRow const* operator->() const;
   bool operator==(const CSVIterator& rhs);
   bool operator!=(const CSVIterator& rhs);

private:
   std::istream* m_str;
   CSVRow m_row;
};

//int main()
//{
//   std::ifstream file("AtomicWeights.csv");
//
//   for (CSVIterator loop(file); loop != CSVIterator(); ++loop) {
//      std::cout << "4th Element(" << (*loop)[3] << ")\n";
//   }
//   file.close();
//
//   return 0;
//}

#endif
