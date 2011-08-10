#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>

using namespace Wt;

class HelloApp : public WApplication
{
		public:
			HelloApp(const WEnvironment &env);
		private:
			WLineEdit *nameEdit_;
			WText *greeting_;

			void greet();
};

HelloApp::HelloApp(const WEnvironment &env) : WApplication(env)
{
		setTitle("sMC webFace by dlg0");

		root()->addWidget(new WText("Enter your name"));

		nameEdit_ = new WLineEdit(root());
		nameEdit_->setFocus();

		WPushButton *b = new WPushButton("Press me",root());
		b->setMargin(5,Left);

		root()->addWidget(new WBreak());
		
		greeting_ = new WText(root());

		b->clicked().connect(this, &HelloApp::greet);
		nameEdit_->enterPressed().connect(this, &HelloApp::greet);
}

void HelloApp::greet()
{
		greeting_->setText("Hello "+ nameEdit_->text());
}

WApplication *createApplication(const WEnvironment &env)
{
		return new HelloApp(env);
}

int main (int argc, char **argv)
{
		return WRun(argc, argv, &createApplication);
}
