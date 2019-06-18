function confidence_table = get_annotation_confidence(filename)
% confidence_table = GET_ANNOTATION_CONFIDENCE(filename)
%
%   filename should be of the form 'run402_annotations.mat'

S = load(filename);
A = S.annotations;
confidence_table = A(A.t==1, {'neuron_id', 'confidence'});
confidence_table.Row = string(confidence_table.neuron_id);