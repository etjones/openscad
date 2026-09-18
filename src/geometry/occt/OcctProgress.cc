#include "geometry/occt/OcctProgress.h"

#include <Message_ProgressIndicator.hxx>
#include <Message_ProgressScope.hxx>
#include <Standard_Handle.hxx>
#include <chrono>

namespace OcctProgress {

namespace {

using Clock = std::chrono::steady_clock;

class Deadline : public Message_ProgressIndicator
{
public:
  explicit Deadline(double seconds)
    : deadline_(Clock::now() +
                std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(seconds)))
  {
  }

  Standard_Boolean UserBreak() override
  {
    if (Clock::now() < deadline_) return Standard_False;
    stopped_ = true;
    return Standard_True;
  }

  [[nodiscard]] bool stopped() const { return stopped_; }

  // Message_ProgressIndicator requires it; the budget shows nothing.
  void Show(const Message_ProgressScope&, const Standard_Boolean) override {}

private:
  Clock::time_point deadline_;
  bool stopped_ = false;
};

// One per export, and only while a budget is set. A Message_ProgressRange
// points into its indicator, so the handle has to outlive every range
// handed out.
thread_local Handle(Deadline) current;

}  // namespace

void setBudget(double seconds)
{
  current = seconds > 0 ? new Deadline(seconds) : Handle(Deadline)();
}

void clear()
{
  current = Handle(Deadline)();
}

Message_ProgressRange range()
{
  return current.IsNull() ? Message_ProgressRange() : current->Start();
}

bool stopped()
{
  return !current.IsNull() && current->stopped();
}

}  // namespace OcctProgress
