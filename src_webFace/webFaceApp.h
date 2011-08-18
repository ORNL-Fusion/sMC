#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>
#include <Wt/WFileUpload>
#include <boost/filesystem.hpp>
#include <Wt/WTable>
#include <Wt/WInteractWidget>
#include <Wt/WSignalMapper>
#include <Wt/WIntValidator>
#include <Wt/WLineEdit>
#include <libconfig.h++>
#include <sstream>
#include <Wt/WDialog>
#include <Wt/WMessageBox>

class webFaceApp : public Wt::WApplication
{
		public:
			webFaceApp(const Wt::WEnvironment &env);
		private:
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_, *spoolFileLocation_;
			Wt::WFileUpload *cfg_upload_;
			Wt::WContainerWidget *logContainer_, *runContainer_, *cfgContainer_;
			Wt::WLineEdit *nRowLineEdit_, *nColLineEdit_, *ZLineEdit_, 
					*amuLineEdit_, *eqdskLineEdit_;
			Wt::WValidator *eqdskValidator_;
			Wt::WPushButton *newRunButton_, *uploadRunButton_, *existingRunButton_;

			void greet();
			void fileUploaded();
			void selectRun(Wt::WText *text);
			void logEntry(const std::string &text);
			void logEntry(const std::string &text, const std::string &color);
			void listRuns();
			void uploadRun();
			void newRun();
};


