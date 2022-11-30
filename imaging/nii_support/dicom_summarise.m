function dicom_summarise( dicomlist )

listentries={};
for n=(1:length(dicomlist))
    listentries{n} = [     ' ' num2str(n) ': '                             ...
                           dicomlist{n}{1}.PatientName.FamilyName ', ' ...
                           dicomlist{n}{1}.PatientName.GivenName       ...
                      ' (' dicomlist{n}{1}.Modality ') ['              ...
                           dicomlist{n}{1}.ProtocolName ' ; '          ...
                           num2str(length(dicomlist{n})) ']'     ];
end;

disp(' ');
disp('DICOM List');
disp('----------');
for n=(1:length(listentries))
    disp(listentries{n});
end;
disp(' ');
