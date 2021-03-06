// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: ReinforcementLearning.proto

#ifndef PROTOBUF_ReinforcementLearning_2eproto__INCLUDED
#define PROTOBUF_ReinforcementLearning_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 2006000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 2006001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)

namespace ReinforcementLearning {

// Internal implementation detail -- do not call these.
void  protobuf_AddDesc_ReinforcementLearning_2eproto();
void protobuf_AssignDesc_ReinforcementLearning_2eproto();
void protobuf_ShutdownFile_ReinforcementLearning_2eproto();

class ReinforcementLearningParameter;
class QLearningSolverParameter;

// ===================================================================

class ReinforcementLearningParameter : public ::google::protobuf::Message {
 public:
  ReinforcementLearningParameter();
  virtual ~ReinforcementLearningParameter();

  ReinforcementLearningParameter(const ReinforcementLearningParameter& from);

  inline ReinforcementLearningParameter& operator=(const ReinforcementLearningParameter& from) {
    CopyFrom(from);
    return *this;
  }

  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }

  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }

  static const ::google::protobuf::Descriptor* descriptor();
  static const ReinforcementLearningParameter& default_instance();

  void Swap(ReinforcementLearningParameter* other);

  // implements Message ----------------------------------------------

  ReinforcementLearningParameter* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const ReinforcementLearningParameter& from);
  void MergeFrom(const ReinforcementLearningParameter& from);
  void Clear();
  bool IsInitialized() const;

  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  ::google::protobuf::Metadata GetMetadata() const;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // optional .ReinforcementLearning.QLearningSolverParameter qLearningSolverParameter = 100;
  inline bool has_qlearningsolverparameter() const;
  inline void clear_qlearningsolverparameter();
  static const int kQLearningSolverParameterFieldNumber = 100;
  inline const ::ReinforcementLearning::QLearningSolverParameter& qlearningsolverparameter() const;
  inline ::ReinforcementLearning::QLearningSolverParameter* mutable_qlearningsolverparameter();
  inline ::ReinforcementLearning::QLearningSolverParameter* release_qlearningsolverparameter();
  inline void set_allocated_qlearningsolverparameter(::ReinforcementLearning::QLearningSolverParameter* qlearningsolverparameter);

  // @@protoc_insertion_point(class_scope:ReinforcementLearning.ReinforcementLearningParameter)
 private:
  inline void set_has_qlearningsolverparameter();
  inline void clear_has_qlearningsolverparameter();

  ::google::protobuf::UnknownFieldSet _unknown_fields_;

  ::google::protobuf::uint32 _has_bits_[1];
  mutable int _cached_size_;
  ::ReinforcementLearning::QLearningSolverParameter* qlearningsolverparameter_;
  friend void  protobuf_AddDesc_ReinforcementLearning_2eproto();
  friend void protobuf_AssignDesc_ReinforcementLearning_2eproto();
  friend void protobuf_ShutdownFile_ReinforcementLearning_2eproto();

  void InitAsDefaultInstance();
  static ReinforcementLearningParameter* default_instance_;
};
// -------------------------------------------------------------------

class QLearningSolverParameter : public ::google::protobuf::Message {
 public:
  QLearningSolverParameter();
  virtual ~QLearningSolverParameter();

  QLearningSolverParameter(const QLearningSolverParameter& from);

  inline QLearningSolverParameter& operator=(const QLearningSolverParameter& from) {
    CopyFrom(from);
    return *this;
  }

  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }

  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }

  static const ::google::protobuf::Descriptor* descriptor();
  static const QLearningSolverParameter& default_instance();

  void Swap(QLearningSolverParameter* other);

  // implements Message ----------------------------------------------

  QLearningSolverParameter* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const QLearningSolverParameter& from);
  void MergeFrom(const QLearningSolverParameter& from);
  void Clear();
  bool IsInitialized() const;

  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  ::google::protobuf::Metadata GetMetadata() const;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // optional int32 numTrainingEpisodes = 1;
  inline bool has_numtrainingepisodes() const;
  inline void clear_numtrainingepisodes();
  static const int kNumTrainingEpisodesFieldNumber = 1;
  inline ::google::protobuf::int32 numtrainingepisodes() const;
  inline void set_numtrainingepisodes(::google::protobuf::int32 value);

  // optional double learningRate = 2 [default = 0.1];
  inline bool has_learningrate() const;
  inline void clear_learningrate();
  static const int kLearningRateFieldNumber = 2;
  inline double learningrate() const;
  inline void set_learningrate(double value);

  // optional double epsilon = 3 [default = 0.95];
  inline bool has_epsilon() const;
  inline void clear_epsilon();
  static const int kEpsilonFieldNumber = 3;
  inline double epsilon() const;
  inline void set_epsilon(double value);

  // optional int32 EpisodeLength = 4;
  inline bool has_episodelength() const;
  inline void clear_episodelength();
  static const int kEpisodeLengthFieldNumber = 4;
  inline ::google::protobuf::int32 episodelength() const;
  inline void set_episodelength(::google::protobuf::int32 value);

  // optional double discount = 5 [default = 0.95];
  inline bool has_discount() const;
  inline void clear_discount();
  static const int kDiscountFieldNumber = 5;
  inline double discount() const;
  inline void set_discount(double value);

  // optional int32 numEpisodesBeforeTraining = 6;
  inline bool has_numepisodesbeforetraining() const;
  inline void clear_numepisodesbeforetraining();
  static const int kNumEpisodesBeforeTrainingFieldNumber = 6;
  inline ::google::protobuf::int32 numepisodesbeforetraining() const;
  inline void set_numepisodesbeforetraining(::google::protobuf::int32 value);

  // optional int32 QTableOutputInterval = 7;
  inline bool has_qtableoutputinterval() const;
  inline void clear_qtableoutputinterval();
  static const int kQTableOutputIntervalFieldNumber = 7;
  inline ::google::protobuf::int32 qtableoutputinterval() const;
  inline void set_qtableoutputinterval(::google::protobuf::int32 value);

  // optional int32 controlInterval = 8 [default = 1];
  inline bool has_controlinterval() const;
  inline void clear_controlinterval();
  static const int kControlIntervalFieldNumber = 8;
  inline ::google::protobuf::int32 controlinterval() const;
  inline void set_controlinterval(::google::protobuf::int32 value);

  // optional int32 experienceReplayInterval = 9 [default = 100];
  inline bool has_experiencereplayinterval() const;
  inline void clear_experiencereplayinterval();
  static const int kExperienceReplayIntervalFieldNumber = 9;
  inline ::google::protobuf::int32 experiencereplayinterval() const;
  inline void set_experiencereplayinterval(::google::protobuf::int32 value);

  // optional int32 experienceStopCriterion = 10 [default = 10000];
  inline bool has_experiencestopcriterion() const;
  inline void clear_experiencestopcriterion();
  static const int kExperienceStopCriterionFieldNumber = 10;
  inline ::google::protobuf::int32 experiencestopcriterion() const;
  inline void set_experiencestopcriterion(::google::protobuf::int32 value);

  // @@protoc_insertion_point(class_scope:ReinforcementLearning.QLearningSolverParameter)
 private:
  inline void set_has_numtrainingepisodes();
  inline void clear_has_numtrainingepisodes();
  inline void set_has_learningrate();
  inline void clear_has_learningrate();
  inline void set_has_epsilon();
  inline void clear_has_epsilon();
  inline void set_has_episodelength();
  inline void clear_has_episodelength();
  inline void set_has_discount();
  inline void clear_has_discount();
  inline void set_has_numepisodesbeforetraining();
  inline void clear_has_numepisodesbeforetraining();
  inline void set_has_qtableoutputinterval();
  inline void clear_has_qtableoutputinterval();
  inline void set_has_controlinterval();
  inline void clear_has_controlinterval();
  inline void set_has_experiencereplayinterval();
  inline void clear_has_experiencereplayinterval();
  inline void set_has_experiencestopcriterion();
  inline void clear_has_experiencestopcriterion();

  ::google::protobuf::UnknownFieldSet _unknown_fields_;

  ::google::protobuf::uint32 _has_bits_[1];
  mutable int _cached_size_;
  double learningrate_;
  ::google::protobuf::int32 numtrainingepisodes_;
  ::google::protobuf::int32 episodelength_;
  double epsilon_;
  double discount_;
  ::google::protobuf::int32 numepisodesbeforetraining_;
  ::google::protobuf::int32 qtableoutputinterval_;
  ::google::protobuf::int32 controlinterval_;
  ::google::protobuf::int32 experiencereplayinterval_;
  ::google::protobuf::int32 experiencestopcriterion_;
  friend void  protobuf_AddDesc_ReinforcementLearning_2eproto();
  friend void protobuf_AssignDesc_ReinforcementLearning_2eproto();
  friend void protobuf_ShutdownFile_ReinforcementLearning_2eproto();

  void InitAsDefaultInstance();
  static QLearningSolverParameter* default_instance_;
};
// ===================================================================


// ===================================================================

// ReinforcementLearningParameter

// optional .ReinforcementLearning.QLearningSolverParameter qLearningSolverParameter = 100;
inline bool ReinforcementLearningParameter::has_qlearningsolverparameter() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void ReinforcementLearningParameter::set_has_qlearningsolverparameter() {
  _has_bits_[0] |= 0x00000001u;
}
inline void ReinforcementLearningParameter::clear_has_qlearningsolverparameter() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void ReinforcementLearningParameter::clear_qlearningsolverparameter() {
  if (qlearningsolverparameter_ != NULL) qlearningsolverparameter_->::ReinforcementLearning::QLearningSolverParameter::Clear();
  clear_has_qlearningsolverparameter();
}
inline const ::ReinforcementLearning::QLearningSolverParameter& ReinforcementLearningParameter::qlearningsolverparameter() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.ReinforcementLearningParameter.qLearningSolverParameter)
  return qlearningsolverparameter_ != NULL ? *qlearningsolverparameter_ : *default_instance_->qlearningsolverparameter_;
}
inline ::ReinforcementLearning::QLearningSolverParameter* ReinforcementLearningParameter::mutable_qlearningsolverparameter() {
  set_has_qlearningsolverparameter();
  if (qlearningsolverparameter_ == NULL) qlearningsolverparameter_ = new ::ReinforcementLearning::QLearningSolverParameter;
  // @@protoc_insertion_point(field_mutable:ReinforcementLearning.ReinforcementLearningParameter.qLearningSolverParameter)
  return qlearningsolverparameter_;
}
inline ::ReinforcementLearning::QLearningSolverParameter* ReinforcementLearningParameter::release_qlearningsolverparameter() {
  clear_has_qlearningsolverparameter();
  ::ReinforcementLearning::QLearningSolverParameter* temp = qlearningsolverparameter_;
  qlearningsolverparameter_ = NULL;
  return temp;
}
inline void ReinforcementLearningParameter::set_allocated_qlearningsolverparameter(::ReinforcementLearning::QLearningSolverParameter* qlearningsolverparameter) {
  delete qlearningsolverparameter_;
  qlearningsolverparameter_ = qlearningsolverparameter;
  if (qlearningsolverparameter) {
    set_has_qlearningsolverparameter();
  } else {
    clear_has_qlearningsolverparameter();
  }
  // @@protoc_insertion_point(field_set_allocated:ReinforcementLearning.ReinforcementLearningParameter.qLearningSolverParameter)
}

// -------------------------------------------------------------------

// QLearningSolverParameter

// optional int32 numTrainingEpisodes = 1;
inline bool QLearningSolverParameter::has_numtrainingepisodes() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void QLearningSolverParameter::set_has_numtrainingepisodes() {
  _has_bits_[0] |= 0x00000001u;
}
inline void QLearningSolverParameter::clear_has_numtrainingepisodes() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void QLearningSolverParameter::clear_numtrainingepisodes() {
  numtrainingepisodes_ = 0;
  clear_has_numtrainingepisodes();
}
inline ::google::protobuf::int32 QLearningSolverParameter::numtrainingepisodes() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.numTrainingEpisodes)
  return numtrainingepisodes_;
}
inline void QLearningSolverParameter::set_numtrainingepisodes(::google::protobuf::int32 value) {
  set_has_numtrainingepisodes();
  numtrainingepisodes_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.numTrainingEpisodes)
}

// optional double learningRate = 2 [default = 0.1];
inline bool QLearningSolverParameter::has_learningrate() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void QLearningSolverParameter::set_has_learningrate() {
  _has_bits_[0] |= 0x00000002u;
}
inline void QLearningSolverParameter::clear_has_learningrate() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void QLearningSolverParameter::clear_learningrate() {
  learningrate_ = 0.1;
  clear_has_learningrate();
}
inline double QLearningSolverParameter::learningrate() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.learningRate)
  return learningrate_;
}
inline void QLearningSolverParameter::set_learningrate(double value) {
  set_has_learningrate();
  learningrate_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.learningRate)
}

// optional double epsilon = 3 [default = 0.95];
inline bool QLearningSolverParameter::has_epsilon() const {
  return (_has_bits_[0] & 0x00000004u) != 0;
}
inline void QLearningSolverParameter::set_has_epsilon() {
  _has_bits_[0] |= 0x00000004u;
}
inline void QLearningSolverParameter::clear_has_epsilon() {
  _has_bits_[0] &= ~0x00000004u;
}
inline void QLearningSolverParameter::clear_epsilon() {
  epsilon_ = 0.95;
  clear_has_epsilon();
}
inline double QLearningSolverParameter::epsilon() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.epsilon)
  return epsilon_;
}
inline void QLearningSolverParameter::set_epsilon(double value) {
  set_has_epsilon();
  epsilon_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.epsilon)
}

// optional int32 EpisodeLength = 4;
inline bool QLearningSolverParameter::has_episodelength() const {
  return (_has_bits_[0] & 0x00000008u) != 0;
}
inline void QLearningSolverParameter::set_has_episodelength() {
  _has_bits_[0] |= 0x00000008u;
}
inline void QLearningSolverParameter::clear_has_episodelength() {
  _has_bits_[0] &= ~0x00000008u;
}
inline void QLearningSolverParameter::clear_episodelength() {
  episodelength_ = 0;
  clear_has_episodelength();
}
inline ::google::protobuf::int32 QLearningSolverParameter::episodelength() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.EpisodeLength)
  return episodelength_;
}
inline void QLearningSolverParameter::set_episodelength(::google::protobuf::int32 value) {
  set_has_episodelength();
  episodelength_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.EpisodeLength)
}

// optional double discount = 5 [default = 0.95];
inline bool QLearningSolverParameter::has_discount() const {
  return (_has_bits_[0] & 0x00000010u) != 0;
}
inline void QLearningSolverParameter::set_has_discount() {
  _has_bits_[0] |= 0x00000010u;
}
inline void QLearningSolverParameter::clear_has_discount() {
  _has_bits_[0] &= ~0x00000010u;
}
inline void QLearningSolverParameter::clear_discount() {
  discount_ = 0.95;
  clear_has_discount();
}
inline double QLearningSolverParameter::discount() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.discount)
  return discount_;
}
inline void QLearningSolverParameter::set_discount(double value) {
  set_has_discount();
  discount_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.discount)
}

// optional int32 numEpisodesBeforeTraining = 6;
inline bool QLearningSolverParameter::has_numepisodesbeforetraining() const {
  return (_has_bits_[0] & 0x00000020u) != 0;
}
inline void QLearningSolverParameter::set_has_numepisodesbeforetraining() {
  _has_bits_[0] |= 0x00000020u;
}
inline void QLearningSolverParameter::clear_has_numepisodesbeforetraining() {
  _has_bits_[0] &= ~0x00000020u;
}
inline void QLearningSolverParameter::clear_numepisodesbeforetraining() {
  numepisodesbeforetraining_ = 0;
  clear_has_numepisodesbeforetraining();
}
inline ::google::protobuf::int32 QLearningSolverParameter::numepisodesbeforetraining() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.numEpisodesBeforeTraining)
  return numepisodesbeforetraining_;
}
inline void QLearningSolverParameter::set_numepisodesbeforetraining(::google::protobuf::int32 value) {
  set_has_numepisodesbeforetraining();
  numepisodesbeforetraining_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.numEpisodesBeforeTraining)
}

// optional int32 QTableOutputInterval = 7;
inline bool QLearningSolverParameter::has_qtableoutputinterval() const {
  return (_has_bits_[0] & 0x00000040u) != 0;
}
inline void QLearningSolverParameter::set_has_qtableoutputinterval() {
  _has_bits_[0] |= 0x00000040u;
}
inline void QLearningSolverParameter::clear_has_qtableoutputinterval() {
  _has_bits_[0] &= ~0x00000040u;
}
inline void QLearningSolverParameter::clear_qtableoutputinterval() {
  qtableoutputinterval_ = 0;
  clear_has_qtableoutputinterval();
}
inline ::google::protobuf::int32 QLearningSolverParameter::qtableoutputinterval() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.QTableOutputInterval)
  return qtableoutputinterval_;
}
inline void QLearningSolverParameter::set_qtableoutputinterval(::google::protobuf::int32 value) {
  set_has_qtableoutputinterval();
  qtableoutputinterval_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.QTableOutputInterval)
}

// optional int32 controlInterval = 8 [default = 1];
inline bool QLearningSolverParameter::has_controlinterval() const {
  return (_has_bits_[0] & 0x00000080u) != 0;
}
inline void QLearningSolverParameter::set_has_controlinterval() {
  _has_bits_[0] |= 0x00000080u;
}
inline void QLearningSolverParameter::clear_has_controlinterval() {
  _has_bits_[0] &= ~0x00000080u;
}
inline void QLearningSolverParameter::clear_controlinterval() {
  controlinterval_ = 1;
  clear_has_controlinterval();
}
inline ::google::protobuf::int32 QLearningSolverParameter::controlinterval() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.controlInterval)
  return controlinterval_;
}
inline void QLearningSolverParameter::set_controlinterval(::google::protobuf::int32 value) {
  set_has_controlinterval();
  controlinterval_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.controlInterval)
}

// optional int32 experienceReplayInterval = 9 [default = 100];
inline bool QLearningSolverParameter::has_experiencereplayinterval() const {
  return (_has_bits_[0] & 0x00000100u) != 0;
}
inline void QLearningSolverParameter::set_has_experiencereplayinterval() {
  _has_bits_[0] |= 0x00000100u;
}
inline void QLearningSolverParameter::clear_has_experiencereplayinterval() {
  _has_bits_[0] &= ~0x00000100u;
}
inline void QLearningSolverParameter::clear_experiencereplayinterval() {
  experiencereplayinterval_ = 100;
  clear_has_experiencereplayinterval();
}
inline ::google::protobuf::int32 QLearningSolverParameter::experiencereplayinterval() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.experienceReplayInterval)
  return experiencereplayinterval_;
}
inline void QLearningSolverParameter::set_experiencereplayinterval(::google::protobuf::int32 value) {
  set_has_experiencereplayinterval();
  experiencereplayinterval_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.experienceReplayInterval)
}

// optional int32 experienceStopCriterion = 10 [default = 10000];
inline bool QLearningSolverParameter::has_experiencestopcriterion() const {
  return (_has_bits_[0] & 0x00000200u) != 0;
}
inline void QLearningSolverParameter::set_has_experiencestopcriterion() {
  _has_bits_[0] |= 0x00000200u;
}
inline void QLearningSolverParameter::clear_has_experiencestopcriterion() {
  _has_bits_[0] &= ~0x00000200u;
}
inline void QLearningSolverParameter::clear_experiencestopcriterion() {
  experiencestopcriterion_ = 10000;
  clear_has_experiencestopcriterion();
}
inline ::google::protobuf::int32 QLearningSolverParameter::experiencestopcriterion() const {
  // @@protoc_insertion_point(field_get:ReinforcementLearning.QLearningSolverParameter.experienceStopCriterion)
  return experiencestopcriterion_;
}
inline void QLearningSolverParameter::set_experiencestopcriterion(::google::protobuf::int32 value) {
  set_has_experiencestopcriterion();
  experiencestopcriterion_ = value;
  // @@protoc_insertion_point(field_set:ReinforcementLearning.QLearningSolverParameter.experienceStopCriterion)
}


// @@protoc_insertion_point(namespace_scope)

}  // namespace ReinforcementLearning

#ifndef SWIG
namespace google {
namespace protobuf {


}  // namespace google
}  // namespace protobuf
#endif  // SWIG

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_ReinforcementLearning_2eproto__INCLUDED
